from contextlib import ExitStack
from pathlib import Path

import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QCheckBox, QDialog, QFileDialog, QHBoxLayout, QLabel, QMessageBox, QPushButton, QVBoxLayout

from .model import coordinates
from .inserts import EMBEDDED_REMARK, LEGACY_EMBEDDED_REMARK


COLORS = {"Upper leaflet backbone": "#79c9a8", "Lower leaflet backbone": "#b09be4", "Polar atoms": "#f5b968",
          "Water (O)": "#62b8ed", "Cations": "#ff8491", "Anions": "#beeb74", "Lipid hydrogens": "#e0e9eb",
          "Embedded molecules": "#f488cc"}


def backbone_segments(residues):
    """Infer display-only heavy-atom bonds strictly within each lipid residue."""
    cache = {}
    segments = {"Upper leaflet backbone": [], "Lower leaflet backbone": []}
    covalent = {"C": .76, "N": .71, "O": .66, "P": 1.07, "S": 1.05}
    for (chain, _, residue), atoms in residues:
        key = (residue, tuple((name, elem) for name, elem, _ in atoms))
        xyz = np.asarray([point for _, _, point in atoms], dtype=np.float32)
        if key not in cache:
            radii = np.asarray([covalent.get(elem, .76) for _, elem, _ in atoms])
            pairs = []
            for i, pos in enumerate(xyz[:-1]):
                distance = np.linalg.norm(xyz[i+1:] - pos, axis=1)
                near = (distance > .55) & (distance <= radii[i] + radii[i+1:] + .30)
                pairs.extend((i, i+1+j) for j in np.flatnonzero(near))
            cache[key] = np.asarray(pairs, dtype=int).reshape(-1, 2)
        segments["Upper leaflet backbone" if chain == "U" else "Lower leaflet backbone"].append(xyz[cache[key]].reshape(-1,3))
    return {k: np.concatenate(parts) for k, parts in segments.items() if parts}


def read_atoms(path):
    groups = {name: [] for name in COLORS}
    residues = []
    bounds = None
    fixed_remaining = 0
    packmol_file = False
    legacy_xy = False
    with Path(path).open(encoding="ascii", errors="replace") as file:
        for line in file:
            if "Built with Packmol" in line:
                packmol_file = True
            if line.startswith((EMBEDDED_REMARK, LEGACY_EMBEDDED_REMARK)):
                fixed_remaining = int(line.split()[-1])
            if line.startswith("CRYST1"):
                dims = [float(line[i:i+9]) for i in (6, 15, 24)]
                bounds = tuple(v for side in dims for v in (0., side))
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            name = line[12:16].strip()
            residue = line[17:21].strip()
            element = line[76:78].strip().upper() or name.lstrip("0123456789")[:1]
            xyz = coordinates(line)
            if packmol_file and bounds and xyz[0] < -0.1 * bounds[1] and xyz[1] < -0.1 * bounds[3]:
                legacy_xy = True
            if fixed_remaining:
                groups["Embedded molecules"].append(xyz)
                fixed_remaining -= 1
                continue
            chain = line[21:22]
            if chain in ("U", "L") and element != "H":
                key = chain, line[22:27], residue
                if not residues or residues[-1][0] != key:
                    residues.append((key, []))
                residues[-1][1].append((name, element, xyz))
            if residue in ("TIP3", "HOH", "WAT", "SOL"):
                if element != "O":
                    continue
                group = "Water (O)"
            elif residue in ("SOD", "NA", "POT", "K"):
                group = "Cations"
            elif residue in ("CLA", "CL"):
                group = "Anions"
            elif element == "H":
                group = "Lipid hydrogens"
            elif element in ("O", "P", "N"):
                group = "Polar atoms"
            else:
                group = "Upper leaflet backbone" if chain == "U" or (chain != "L" and xyz[2] >= 0) else "Lower leaflet backbone"
            groups[group].append(xyz)
    if not any(groups.values()):
        raise ValueError("The PDB contains no atoms to display.")
    # Standalone centered Packmol files predate corner-based export. The
    # original system.json takes precedence when the exact bounds are known.
    if legacy_xy:
        bounds = tuple(v for side in bounds[1::2] for v in (-side/2, side/2))
    return groups, bounds, backbone_segments(residues)


def save_cell_aligned_copy(source, destination, bounds):
    """Translate a legacy builder PDB using the exact lower corner from its manifest."""
    corner = bounds[::2]
    destination = Path(destination)
    with Path(source).open("rb") as original, ExitStack() as rollback:
        with destination.open("xb") as result:
            rollback.callback(destination.unlink)
            for line in original:
                if line.startswith((b"ATOM  ", b"HETATM")):
                    translated = b"".join(
                        f"{float(line[i:i+8]) - low:8.3f}".encode("ascii")
                        for i, low in zip((30, 38, 46), corner)
                    )
                    if len(translated) != 24:
                        raise ValueError("A translated coordinate exceeds the PDB field width.")
                    line = line[:30] + translated + line[54:]
                result.write(line)
        rollback.pop_all()


class Viewer(QDialog):
    def __init__(self, parent, plan=None, path=None, bounds=None):
        super().__init__(parent)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.setWindowTitle("Packing Constraints in 3D" if plan else f"Packed Bilayer in 3D  ·  {Path(path).name}")
        self.resize(870, 640)
        root = QVBoxLayout(self)
        self.plot = QtInteractor(self)
        root.addWidget(self.plot.interactor, 1)
        self.plot.set_background("#182731")
        controls = QHBoxLayout()
        root.addLayout(controls)
        legacy_bounds = bounds if bounds and any(bounds[::2]) else None
        if plan:
            bounds = plan.cell_bounds
            x0, x1, y0, y1, low, high = bounds
            midplane = plan.membrane_origin[2]
            for start, end in [(midplane + plan.interfaces[0], high), (low, midplane - plan.interfaces[1])]:
                if end > start:
                    self.plot.add_mesh(pv.Box(bounds=(x0,x1,y0,y1,start,end)), color="#62b8ed", opacity=0.22)
            for (_, u, l), h in zip(plan.rows, plan.heads):
                for count, z in [(u, h), (l, -h)]:
                    if count:
                        self.plane(plan, z, "#f5b968")
            for z in (plan.settings.tail_band, -plan.settings.tail_band):
                self.plane(plan, z, "#79c9a8")
            legend = "Blue — water layers     Orange — head boundaries     Green — tail boundaries"
            for placement in plan.placements:
                self.plot.add_points(placement.xyz + plan.membrane_origin,color=COLORS["Embedded molecules"],point_size=7,render_points_as_spheres=True)
            if plan.placements:
                legend += "     Pink — embedded molecules"
            controls.addWidget(QLabel(legend))
        else:
            groups, file_bounds, backbones = read_atoms(path)
            bounds = bounds or file_bounds
            for name, points in groups.items():
                if not points:
                    continue
                size = 10 if name in ("Cations", "Anions") else (3 if name in ("Water (O)", "Lipid hydrogens") else 5)
                actor = self.plot.add_points(np.asarray(points, dtype=np.float32), color=COLORS[name],
                                             point_size=size, render_points_as_spheres=True, render=False)
                actors = [actor]
                if name in backbones and len(backbones[name]):
                    actors.append(self.plot.add_lines(backbones[name], color=COLORS[name], width=2))
                check = QCheckBox(name)
                active = name != "Lipid hydrogens"
                check.setChecked(active)
                check.setStyleSheet(f"color: {COLORS[name]}; background: #24333c; padding: 3px")
                actor.SetVisibility(active)
                check.toggled.connect(lambda on, items=actors: self.visibility(items, on))
                controls.addWidget(check)
        if bounds:
            self.plot.add_mesh(pv.Box(bounds=bounds), style="wireframe", color="#b5c6cf", line_width=1)
        self.plot.add_axes()
        self.plot.camera_position = [(1.5,-2.4,.95), (0,0,0), (0,0,1)]
        self.plot.reset_camera()
        tools = QHBoxLayout()
        root.addLayout(tools)
        tools.addWidget(QLabel("Drag with the left mouse button to rotate. Scroll to zoom."))
        tools.addStretch()
        if path and legacy_bounds:
            save = QPushButton("Save Cell-Aligned PDB…")
            save.clicked.connect(lambda: self.save_aligned(path, legacy_bounds))
            tools.addWidget(save)
        for text, slot in [("Side View", self.plot.view_xz), ("Top View", self.plot.view_xy)]:
            btn = QPushButton(text)
            btn.clicked.connect(lambda _, callback=slot: callback())
            tools.addWidget(btn)

    def plane(self, plan, z, color):
        x, y, midplane = plan.membrane_origin
        mesh = pv.Plane(center=(x,y,midplane + z), direction=(0,0,1), i_size=plan.x, j_size=plan.y,
                        i_resolution=1, j_resolution=1)
        self.plot.add_mesh(mesh, color=color, opacity=.34, show_edges=True, edge_color=color)

    def save_aligned(self, path, bounds):
        target, _ = QFileDialog.getSaveFileName(self, "Save a cell-aligned copy",
            str(Path(path).with_name("bilayer-cell-aligned.pdb")), "PDB (*.pdb)")
        if target:
            try:
                save_cell_aligned_copy(path, target, bounds)
            except FileExistsError:
                QMessageBox.warning(self, "File already exists", "Choose a new filename for the cell-aligned copy.")
            except (OSError, ValueError) as exc:
                QMessageBox.warning(self, "Could not save the PDB", str(exc))

    def visibility(self, actors, on):
        for actor in actors:
            actor.SetVisibility(on)
        self.plot.render()

    def closeEvent(self, event):
        self.plot.close()
        super().closeEvent(event)
