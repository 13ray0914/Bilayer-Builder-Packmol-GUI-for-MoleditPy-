"""Bilayer geometry in Å. Packmol atom indices are one-based file positions."""
from contextlib import ExitStack
from dataclasses import dataclass
from decimal import Decimal, ROUND_HALF_UP
from itertools import count
from math import acos, dist, floor, isfinite, pi, sqrt
from pathlib import Path, PureWindowsPath
import json
import ntpath
import random
import tempfile


TEMPLATES = Path(__file__).parent / "templates"
AVOGADRO_A3 = 6.02214076e-4  # molecules in one Å³ at one mol/L
WATER_MOLARITY = 1000 / 18.01528  # pure water at 1 g/cm³

LIPID_GROUPS = {
    "Sterols": "Sterols",
    "PA": "PA (Phosphatidic Acid) Lipids",
    "PC": "PC (Phosphatidylcholine) Lipids",
    "PE": "PE (Phosphatidylethanolamine) Lipids",
    "PG": "PG (Phosphatidylglycerol) Lipids",
    "PS": "PS (Phosphatidylserine) Lipids",
    "Custom": "Custom Lipids",
}

# Acyl chains in sn-1 / sn-2 order for the bundled phospholipid templates.
CATALOG_TAILS = {
    "DL": "12:0 / 12:0", "DM": "14:0 / 14:0", "DP": "16:0 / 16:0",
    "DS": "18:0 / 18:0", "PO": "16:0 / 18:1", "DO": "18:1 / 18:1",
    "SD": "18:0 / 22:6", "DA": "20:4 / 20:4",
}


def coordinates(line):
    xyz = tuple(float(line[i:i + 8]) for i in (30, 38, 46))
    if not all(map(isfinite, xyz)):
        raise ValueError("The PDB contains a non-finite coordinate.")
    return xyz


def molecule_records(path):
    lines = Path(path).read_text(encoding="ascii", errors="replace").splitlines()
    atoms = [s for s in lines if s.startswith(("ATOM  ", "HETATM"))]
    if not atoms or len({s[17:27] for s in atoms}) != 1 or sum(s.startswith("MODEL ") for s in lines) > 1:
        raise ValueError("Choose a PDB containing one molecule, one residue and one model.")
    if any(s[16:17].strip() for s in atoms):
        raise ValueError("Keep one conformation and remove alternate-location indicators from the PDB.")
    for line in atoms:
        coordinates(line)
    return tuple(atoms)


@dataclass(frozen=True)
class Lipid:
    name: str
    atoms: tuple
    head: int
    tails: tuple
    area: float
    charge: int
    group: str = "Custom"
    tail_info: str = "Unspecified"

    def __post_init__(self):
        marks = (self.head,) + self.tails
        if not self.tails or len(set(marks)) != len(marks) or not all(1 <= i <= len(self.atoms) for i in marks):
            raise ValueError("Select distinct, valid atom indices for the head and each tail end.")
        if self.length <= 1 or self.area <= 0:
            raise ValueError("The head-to-tail distance must exceed 1 Å and the area per lipid must be positive.")

    @property
    def length(self):
        ends = [coordinates(self.atoms[i - 1]) for i in self.tails]
        mean_end = tuple(sum(axis) / len(ends) for axis in zip(*ends))
        return dist(coordinates(self.atoms[self.head - 1]), mean_end)


def library():
    def catalog_lipid(name, path, head, tails, area, charge):
        if name == "CHL1":
            group, tail_info = "Sterols", "N/A"
        else:
            group, tail_info = name[2:], CATALOG_TAILS[name[:2]]
        return Lipid(name, molecule_records(path), head, tails, area, charge, group, tail_info)

    specs = [
        ("POPC", (88, 131), 64, 0), ("POPE", (79, 122), 57, 0),
        ("POPS", (81, 124), 63, -1), ("POPG", (81, 124), 66, -1),
        ("DPPC", (84, 127), 63, 0), ("DOPC", (88, 135), 67, 0),
        ("DMPC", (78, 115), 61, 0),
    ]
    lipids = [catalog_lipid(name, TEMPLATES / (name.lower() + ".pdb"), 1, tails, area, charge)
              for name, tails, area, charge in specs]
    extra = TEMPLATES / "memgen"
    for item in json.loads((extra / "lipids.json").read_text(encoding="ascii")):
        name = item["name"]
        lipids.append(catalog_lipid(name, extra / (name + ".pdb"), item["head"],
                            tuple(item["tails"]), item["area"], item["charge"]))
    return lipids


@dataclass
class Settings:
    water: float = 18.0  # water slab thickness on each side
    salt: float = 0.15
    cation: str = "Na"
    neutralize: bool = True
    xy: tuple = ()  # empty selects composition-based square
    area_scale: float = 1.20  # initial packing area relative to the reference areas
    extension: float = 0.88  # fraction of head-to-tail distance projected onto z
    tail_band: float = 4.0
    tolerance: float = 2.0
    hydrogen_scale: float = 0.7
    seed: int = 12345
    loops: int = 200
    maxit: int = 200
    lattice: bool = True
    water_mode: str = "volume"
    waters_per_lipid: float = 30.0
    total_waters: int = 1000


def create_output_folder(parent, name):
    if not name.strip():
        return Path(tempfile.mkdtemp(prefix="bilayer-", dir=parent))
    if (name.endswith((" ", ".")) or any(c in '<>:"/\\|?*' or ord(c) < 32 for c in name)
            or (ntpath.isreserved(name) if hasattr(ntpath, "isreserved") else PureWindowsPath(name).is_reserved())):
        raise ValueError("Enter a single valid output folder name, such as POPC-run-01.")
    for index in count(1):
        folder = Path(parent) / (name if index == 1 else f"{name}-{index}")
        try:
            folder.mkdir()
        except FileExistsError:
            continue
        return folder


def layer_grid(count, width, depth):
    if not count:
        return []
    nrows = min(count, max(1, round(sqrt(count * depth / width))))
    points = []
    for row in range(nrows):
        ncols = count // nrows + (row < count % nrows)
        for col in range(ncols):
            points.append((width*((col+.5)/ncols-.5), depth*((row+.5)/nrows-.5)))
    return points


class Plan:
    def __init__(self, rows, settings, placements=()):
        self.settings = s = settings
        self.placements = tuple(placements)
        self.excluded = None
        self.rows = [(mol, u, l) for mol, u, l in rows if u or l]
        self.counts = [sum(row[i] for row in self.rows) for i in (1, 2)]
        if min(self.counts) < 1:
            raise ValueError("Add at least one lipid molecule to each leaflet.")
        self.areas = [sum(row[0].area * row[i] for row in self.rows) for i in (1, 2)]
        self.heads = [s.tail_band / 2 + mol.length * s.extension for mol, _, _ in self.rows]
        # Each leaflet can have a different longest species.
        self.interfaces = [max(h for h, row in zip(self.heads, self.rows) if row[i]) + 2 for i in (1, 2)]
        side = sqrt(max(self.areas) * s.area_scale)
        if self.placements:
            from .inserts import ExcludedSpace
            self.excluded = excluded = ExcludedSpace(self.placements, self.interfaces, s.tolerance)
            side = max(sqrt(area*s.area_scale + len(mask)*4) for area,mask in zip(self.areas,excluded.footprints))
            side = max(side, 2*max(excluded.extents) + s.tolerance + .01)
        self.x, self.y = s.xy or (side, side)
        self.bounds = (-self.x/2, self.x/2, -self.y/2, self.y/2,
                       -self.interfaces[1] - s.water, self.interfaces[0] + s.water)
        volumes = [self.x * self.y * s.water] * 2
        if self.excluded:
            if any(2*r + s.tolerance > length for r,length in zip(excluded.extents,(self.x,self.y))):
                raise ValueError("An embedded structure crosses the XY cell boundary. Enlarge the cell or reposition it.")
            if excluded.zmin < self.bounds[4]+s.tolerance/2 or excluded.zmax > self.bounds[5]-s.tolerance/2:
                raise ValueError("An embedded structure crosses the Z boundary. Increase water thickness or reposition it.")
            volumes = excluded.solvent_volumes(self.x*self.y, s.water, self.interfaces)
        self.volume = sum(volumes)
        self.pairs = floor(s.salt * self.volume * AVOGADRO_A3 + 0.5)
        self.charge = sum(mol.charge * (u + l) for mol, u, l in self.rows) + sum(p.charge for p in self.placements)
        extra = self.charge if s.neutralize else 0
        self.cations = self.pairs + max(0, -extra)
        self.anions = self.pairs + max(0, extra)
        # Place salt in pairs. Split neutralizing ions in proportion to the
        # like-signed lipid charge of the two leaflets.
        upper_fraction = volumes[0]/self.volume if self.volume else .5
        upper_pairs = floor(self.pairs*upper_fraction+.5)
        leaflet_charge = [sum(row[0].charge * row[i] for row in self.rows) for i in (1, 2)]
        for placement in self.placements:
            fraction = (placement.xyz[:,2] >= 0).mean()
            leaflet_charge[0] += placement.charge*fraction
            leaflet_charge[1] += placement.charge*(1-fraction)
        weights = [max(0, q if extra > 0 else -q) for q in leaflet_charge]
        upper_extra = floor(abs(extra) * weights[0] / sum(weights) + .5) if extra else 0
        upper_cations = upper_pairs + (upper_extra if extra < 0 else 0)
        upper_anions = upper_pairs + (upper_extra if extra > 0 else 0)
        if s.water_mode == "volume":
            sites = floor(self.volume * AVOGADRO_A3 * WATER_MOLARITY + .5)
            self.water_count = sites - self.cations - self.anions
            upper_water = floor(sites * upper_fraction + .5) - upper_cations - upper_anions
        else:
            if s.water_mode == "per_lipid":
                if not isfinite(s.waters_per_lipid) or s.waters_per_lipid < 0:
                    raise ValueError("Water molecules per lipid must be a finite, nonnegative number.")
                target = Decimal(str(s.waters_per_lipid)) * sum(self.counts)
                self.water_count = int(target.to_integral_value(rounding=ROUND_HALF_UP))
            elif s.water_mode == "total":
                if not isinstance(s.total_waters, int) or s.total_waters < 0:
                    raise ValueError("Total water molecules must be a nonnegative integer.")
                self.water_count = s.total_waters
            else:
                raise ValueError("Choose a water count method.")
            if self.water_count and self.volume <= 0:
                raise ValueError("Increase the water thickness to accommodate the requested water molecules.")
            if self.water_count == 0 and (s.salt or self.cations or self.anions):
                raise ValueError("For zero water molecules, set added salt to zero and turn off any required counterions.")
            upper_water = floor(self.water_count * upper_fraction + .5)
        self.slabs = {"water": (upper_water, self.water_count - upper_water),
                      "cation": (upper_cations, self.cations - upper_cations),
                      "chloride": (upper_anions, self.anions - upper_anions)}
        if min(self.slabs["water"]) < 0 or (s.water == 0 and (self.cations or self.anions or s.salt)):
            raise ValueError("Increase the water thickness to accommodate the salt or counterions.")
        self.atom_count = sum(len(m.atoms)*(u+l) for m, u, l in self.rows) + 3*self.water_count + self.cations + self.anions
        self.atom_count += sum(len(p.xyz) for p in self.placements)

    @property
    def summary(self):
        s = self.settings
        conc = self.pairs / (self.volume * AVOGADRO_A3) if self.volume else 0
        z = self.bounds[5] - self.bounds[4]
        embedded = f"     ·     Embedded structures {len(self.placements)}" if self.placements else ""
        return (f"Lipids per leaflet: {self.counts[0]} upper / {self.counts[1]} lower     ·     "
                f"Cell {self.x:.1f} × {self.y:.1f} × {z:.1f} Å{embedded}\n"
                f"Water {self.water_count:,} ({self.slabs['water'][0]:,} upper / {self.slabs['water'][1]:,} lower)     "
                f"{s.cation}⁺ {self.cations}     Cl⁻ {self.anions}     "
                f"Added salt {conc:.3f} M     Ion pairs {self.pairs}     Net charge {self.charge+self.cations-self.anions:+d} e")

    def deck(self):
        s = self.settings
        a, b, c, d, low, high = self.bounds
        lines = ["filetype pdb", "output bilayer.pdb", "packall", "movebadrandom", f"maxit {s.maxit}", f"tolerance {s.tolerance:.3f}",
                 f"seed {s.seed}", f"nloop {s.loops}", f"pbc {a:.4f} {c:.4f} {low:.4f} {b:.4f} {d:.4f} {high:.4f}"]
        if s.lattice:
            lines.append("restart_from placement.dat")
        if self.placements:
            lines.append("amber_ter_preserve")
        for i in range(1,len(self.placements)+1):
            lines.extend(["", f"structure embedded_{i}.pdb", "  number 1", "  resnumbers 1",
                          "  fixed 0. 0. 0. 0. 0. 0.", "end structure"])
        def hydrogen_radius(indices):
            if s.hydrogen_scale == 1:
                return []
            result = []
            for start in range(0, len(indices), 12):
                group = indices[start:start+12]
                result.extend(["  atoms " + " ".join(map(str, group)),
                               f"    radius {s.tolerance * s.hydrogen_scale / 2:.4f}", "  end atoms"])
            return result

        def block(file, count, chain, z0, z1, constraints=()):
            if count:
                lines.extend(["", f"structure {file}.pdb", f"  number {count}", f"  chain {chain}", "  resnumbers 3",
                              f"  inside box {a:.4f} {c:.4f} {z0:.4f} {b:.4f} {d:.4f} {z1:.4f}",
                              *constraints, "end structure"])
        for i, ((mol, upper, lower), head) in enumerate(zip(self.rows, self.heads)):
            for sign, count, chain, interface in [(1, upper, "U", self.interfaces[0]), (-1, lower, "L", self.interfaces[1])]:
                outer = "over" if sign == 1 else "below"
                inner = "below" if sign == 1 else "over"
                planes = [f"  atoms {mol.head}", f"    {outer} plane 0 0 1 {sign * head:.4f}", "  end atoms",
                          "  atoms " + " ".join(map(str, mol.tails)),
                          f"    {inner} plane 0 0 1 {sign * s.tail_band:.4f}", "  end atoms"]
                planes.extend(hydrogen_radius([j for j, atom in enumerate(mol.atoms, 1) if atom[76:78].strip() == "H"]))
                z0, z1 = sorted((-sign * s.tail_band, sign * (interface + 3)))
                block(f"lipid_{i+1}", count, chain, max(low, z0), min(high, z1), planes)
        for file, (upper, lower) in self.slabs.items():
            radii = hydrogen_radius([2, 3]) if file == "water" else []
            block(file, upper, "A", self.interfaces[0], high, radii)
            block(file, lower, "B", low, -self.interfaces[1], radii)
        return "\n".join(lines) + "\n"

    def write_placement(self, stream):
        # Packmol restart entries are initial translation and three Euler angles
        # in radians, in structure-block order. All six variables remain free.
        rng = random.Random(self.settings.seed)
        def write(*values):
            stream.write(" ".join(f"{v:.7f}" for v in values) + "\n")

        def free_layer(n, side, z=None):
            points = layer_grid(n, self.x, self.y)
            if not self.excluded or n == 0:
                return points
            count = n
            for _ in range(12):
                free = [p for p in points if not self.excluded.blocked(*p,side,z)]
                if len(free) >= n:
                    return rng.sample(free,n)
                count += max(1,n-len(free),count//10)
                points = layer_grid(count,self.x,self.y)
            raise ValueError("Too few free grid sites around the embedded structure. Enlarge the cell or disable grid placement.")

        grids = [free_layer(n,side) for side,n in enumerate(self.counts)]
        for grid in grids:
            rng.shuffle(grid)
        from .seeding import LipidSeedGrid
        seed_grid = LipidSeedGrid(self)
        for (mol, u, l), boundary in zip(self.rows, self.heads):
            for side, count in enumerate((u,l)):
                candidates = seed_grid.candidates(mol,side,boundary,rng)
                for _ in range(count):
                    x, y = grids[side].pop()
                    write(*seed_grid.choose(x,y,candidates,rng))
        solvent = [[], []]
        for side, sign in enumerate((1,-1)):
            total = sum(pair[side] for pair in self.slabs.values())
            if not total:
                continue
            nz = min(total, max(1, round(self.settings.water / 3.10)))
            for iz in range(nz):
                n = total//nz + (iz < total % nz)
                z = sign*(self.interfaces[side] + self.settings.water*(iz+.5)/nz)
                solvent[side].extend((x,y,z) for x,y in free_layer(n,side,z))
            rng.shuffle(solvent[side])
        for counts in self.slabs.values():
            for side, count in enumerate(counts):
                for _ in range(count):
                    write(*solvent[side].pop(), rng.uniform(0,2*pi), rng.uniform(0,2*pi), acos(rng.uniform(-1,1)))

    def export(self, parent, executable=None, *, folder_name=""):
        folder = create_output_folder(parent, folder_name)
        with ExitStack() as rollback:
            rollback.callback(folder.rmdir)
            def output_file(name):
                path = folder / name
                rollback.callback(path.unlink, missing_ok=True)
                return path
            self._write_inputs(output_file, executable)
            rollback.pop_all()
        return folder

    def _write_inputs(self, output_file, executable):
        def write_template(name, atoms):
            output_file(name).write_text("\n".join(atoms) + "\nEND\n", encoding="ascii")
        manifest = []
        embedded = []
        for i,p in enumerate(self.placements,1):
            name = f"embedded_{i}.pdb"
            p.write(output_file(name))
            embedded.append(dict(name=p.structure.name, pdb=name, atoms=len(p.xyz), charge=p.charge,
                                 shift=p.shift, rotation_xyz_deg=p.angles, center=p.center))
        for i, (mol, u, l) in enumerate(self.rows, 1):
            name = f"lipid_{i}.pdb"
            write_template(name, mol.atoms)
            manifest.append(dict(name=mol.name, pdb=name, upper=u, lower=l, head=mol.head, tails=mol.tails,
                                 area=mol.area, charge=mol.charge, group=mol.group, tail_info=mol.tail_info))
        for name in ("water", "chloride", "sodium"):
            if name == "sodium" and self.settings.cation == "K":
                write_template("cation.pdb", ["HETATM    1  POT POT X   1       0.000   0.000   0.000  1.00  0.00           K"])
            else:
                target = "cation" if name == "sodium" else name
                write_template(target + ".pdb", molecule_records(TEMPLATES / (name + ".pdb")))
        output_file("packmol.inp").write_text(self.deck(), encoding="ascii")
        if self.settings.lattice:
            with output_file("placement.dat").open("w", encoding="ascii") as stream:
                self.write_placement(stream)
        output_file("system.json").write_text(json.dumps(dict(settings=vars(self.settings), lipids=manifest,
            executable=str(executable) if executable else None, embedded=embedded,
            bounds=self.bounds, water=self.water_count, water_by_leaflet=self.slabs["water"],
            cations=self.cations, anions=self.anions,
            added_pairs=self.pairs, atoms=self.atom_count), ensure_ascii=False, indent=2), encoding="utf-8")
