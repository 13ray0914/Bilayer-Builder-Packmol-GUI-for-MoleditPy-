"""Rigid PDB/SDF inclusions and geometry estimates for initial packing."""
from dataclasses import dataclass
from functools import cached_property
from itertools import product
from math import floor
from pathlib import Path

import numpy as np
import os
import shutil
import tempfile


from . import PLUGIN_NAME


EMBEDDED_REMARK = f"REMARK 950 {PLUGIN_NAME} FIXED ATOMS "
LEGACY_EMBEDDED_REMARK = "REMARK 950 MD BILAYER BUILDER FIXED ATOMS "


def verify_fixed_output(path, placements):
    """Check the fixed prefix and embed its length so the PDB travels on its own."""
    if not placements:
        return
    n = sum(len(p.xyz) for p in placements)
    pending = None
    try:
        with path.open("rb") as source, tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as out:
            pending = Path(out.name)
            out.write(f"{EMBEDDED_REMARK}{n}\n".encode("ascii"))
            for p in placements:
                for original, xyz in zip(p.structure.atoms,p.xyz):
                    for raw in source:
                        out.write(raw)
                        if raw.startswith((b"ATOM  ",b"HETATM")):
                            break
                    else:
                        raise ValueError("Packmol output is missing embedded atoms.")
                    line = raw.decode("ascii")
                    packed = [float(line[i:i+8]) for i in (30,38,46)]
                    if (line[12:27] != original[12:27] or
                            np.max(np.abs(packed-xyz)) > .002):
                        raise ValueError(f"Packmol changed the placement or atom identity of {p.structure.name}.")
            shutil.copyfileobj(source,out)
        os.replace(pending,path)
    finally:
        if pending is not None:
            pending.unlink(missing_ok=True)


@dataclass
class Structure:
    name: str
    atoms: list
    coords: np.ndarray
    connections: list
    breaks: set
    charge: int

    @classmethod
    def load(cls, path):
        path = Path(path)
        charge = 0
        if path.suffix.lower() in (".sdf", ".mol"):
            from rdkit import Chem
            if path.suffix.lower() == ".sdf":
                reader = Chem.SDMolSupplier(str(path), removeHs=False)
                if len(reader) != 1:
                    raise ValueError("Select an SDF containing exactly one molecule.")
                mol = reader[0]
            else:
                mol = Chem.MolFromMolFile(str(path), removeHs=False)
            if mol is None or not mol.GetNumConformers() or not mol.GetConformer().Is3D():
                raise ValueError("Supply a MOL or SDF with a 3D conformer.")
            charge = Chem.GetFormalCharge(mol)
            lines = Chem.MolToPDBBlock(mol).splitlines()
        else:
            lines = path.read_text(encoding="ascii", errors="replace").splitlines()
            if sum(line.startswith("MODEL ") for line in lines) > 1:
                raise ValueError("Select one PDB model before embedding.")
        atoms, coords, connections, breaks, serials = [], [], [], set(), {}
        for line in lines:
            if line.startswith(("ATOM  ", "HETATM")):
                if line[16:17].strip():
                    raise ValueError("Resolve alternate locations before embedding.")
                serial = int(line[6:11])
                if serial in serials:
                    raise ValueError("Embedded PDB atom serials must be unique.")
                int(line[22:26])  # Packmol requires numeric residue numbers.
                coords.append([float(line[i:i+8]) for i in (30,38,46)])
                atoms.append(line)
                serials[serial] = len(atoms)
            elif line.startswith("TER"):
                breaks.add(len(atoms))
        if not atoms or not np.isfinite(coords).all():
            raise ValueError("The embedded molecule must contain atoms with finite 3D coordinates.")
        for line in lines:
            if line.startswith("CONECT"):
                fields = [line[i:i+5] for i in range(6,len(line.rstrip()),5)]
                try:
                    ids = [serials[int(f)] for f in fields if f.strip()]
                except KeyError as exc:
                    raise ValueError("A CONECT record refers to a missing atom.") from exc
                if len(ids) > 1:
                    connections.append(ids)
        return cls(path.name, atoms, np.array(coords), connections, breaks, charge)


@dataclass
class Placement:
    structure: Structure
    shift: tuple
    angles: tuple
    center: bool
    charge: int

    @cached_property
    def xyz(self):
        xyz = self.structure.coords.copy()
        if self.center:
            xyz -= xyz.mean(axis=0)
        for angle, (a,b) in zip(np.radians(self.angles), ((1,2),(2,0),(0,1))):
            c, s = np.cos(angle), np.sin(angle)
            xyz[:, [a,b]] = xyz[:, [a,b]] @ np.array([[c,s],[-s,c]])
        return xyz + np.asarray(self.shift)

    def write(self, path):
        with path.open("w", encoding="ascii") as out:
            for i,(line,xyz) in enumerate(zip(self.structure.atoms,self.xyz),1):
                out.write(f"{line[:6]}{i:5d}{line[11:30]}" + "".join(f"{v:8.3f}" for v in xyz) + line[54:] + "\n")
                if i in self.structure.breaks:
                    out.write("TER\n")
            for ids in self.structure.connections:
                out.write("CONECT" + "".join(f"{i:5d}" for i in ids) + "\n")
            out.write("END\n")


class ExcludedSpace:
    """Estimate occupied area/volume with 2 Å voxels around fixed atomic centers."""
    def __init__(self, placements, interfaces, tolerance):
        self.cells = set()
        self.footprints = [set(), set()]
        all_points = np.concatenate([p.xyz for p in placements])
        self.extents = np.max(np.abs(all_points[:,:2]), axis=0)
        self.zmin, self.zmax = all_points[:,2].min(), all_points[:,2].max()
        # Check fixed/fixed distances; Packmol cannot resolve these overlaps.
        prior = {}
        neighbors = tuple(product((-1,0,1), repeat=3))
        for p in placements:
            xyz = p.xyz
            buckets = np.floor(xyz / tolerance).astype(int)
            if prior:
                for atom, (x,y,z) in zip(xyz,buckets):
                    for dx,dy,dz in neighbors:
                        other = prior.get((x+dx,y+dy,z+dz))
                        if other is not None and np.any(np.sum((other-atom)**2,axis=1) < tolerance**2):
                            raise ValueError("Fixed molecules overlap. Adjust their positions before building.")
            if len(placements) > 1:
                for atom, idx in zip(xyz,buckets):
                    prior.setdefault(tuple(idx), []).append(atom)
        base = np.floor(all_points / 2).astype(int)
        for delta in neighbors:
            indices = base + delta
            within = np.sum(((indices+.5)*2 - all_points)**2, axis=1) <= 4.
            self.cells.update(map(tuple,indices[within].tolist()))
        for x,y,z in self.cells:
            height = 2*z + 1
            if 0 <= height <= interfaces[0]:
                self.footprints[0].add((x,y))
            elif -interfaces[1] <= height < 0:
                self.footprints[1].add((x,y))

    def solvent_volumes(self, area, water, interfaces):
        volumes = [area*water, area*water]
        for _,_,k in self.cells:
            for side,(lo,hi) in enumerate(((interfaces[0],interfaces[0]+water),(-interfaces[1]-water,-interfaces[1]))):
                volumes[side] -= 4*max(0.,min(2*k+2,hi)-max(2*k,lo))
        return [max(0.,v) for v in volumes]

    def blocked(self, x, y, side, z):
        xy = floor(x/2), floor(y/2)
        return xy in self.footprints[side] if z is None else (*xy,floor(z/2)) in self.cells
