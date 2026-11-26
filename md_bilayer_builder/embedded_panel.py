from PySide6.QtCore import Signal
from PySide6.QtWidgets import (QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QFileDialog, QFormLayout, QHBoxLayout, QLabel, QListWidget, QMessageBox,
    QPushButton, QSpinBox, QVBoxLayout, QWidget)

from .inserts import Placement, Structure


class EmbeddedPanel(QWidget):
    changed = Signal()

    def __init__(self):
        super().__init__()
        self.placements = []
        layout = QVBoxLayout(self)
        hint = QLabel("Add a membrane protein, complex or synthetic molecule. Packmol holds each structure "
                      "at its chosen position while packing lipids, water and ions around it.")
        hint.setWordWrap(True)
        layout.addWidget(hint)
        self.list = QListWidget()
        layout.addWidget(self.list)
        bar = QHBoxLayout()
        layout.addLayout(bar)
        for label,slot in [("Add Structure…",self.add), ("Edit Placement…",self.edit), ("Remove",self.remove)]:
            btn = QPushButton(label)
            btn.clicked.connect(slot)
            bar.addWidget(btn)
        bar.addStretch()
        self.list.itemDoubleClicked.connect(lambda *_: self.edit())

    def add(self):
        path,_ = QFileDialog.getOpenFileName(self,"Open a structure to embed","","Structures (*.pdb *.sdf *.mol);;All files (*)")
        if path:
            try:
                structure = Structure.load(path)
                self.set_placement(Placement(structure,(0.,)*3,(0.,)*3,True,structure.charge))
            except (ValueError,OSError) as exc:
                QMessageBox.warning(self,"Embedded Molecules",str(exc))

    def edit(self):
        row = self.list.currentRow()
        if row >= 0:
            self.set_placement(self.placements[row],row)

    def remove(self):
        row = self.list.currentRow()
        if row >= 0:
            self.placements.pop(row)
            self.list.takeItem(row)
            self.changed.emit()

    def set_placement(self,p,row=None):
        dlg = QDialog(self)
        dlg.setWindowTitle(f"Placement  ·  {p.structure.name}")
        form = QFormLayout(dlg)
        center = QCheckBox("Move the geometric center to the origin first")
        center.setChecked(p.center)
        form.addRow(center)
        fields = []
        for label,values,limit in (("Rotate (°)",p.angles,360),("Translate (Å)",p.shift,2000)):
            group = QHBoxLayout()
            axes = []
            for axis,value in zip("XYZ",values):
                spin = QDoubleSpinBox()
                spin.setRange(-limit,limit)
                spin.setDecimals(2)
                spin.setValue(value)
                group.addWidget(QLabel(axis))
                group.addWidget(spin)
                axes.append(spin)
            form.addRow(label,group)
            fields.append(axes)
        charge = QSpinBox()
        charge.setRange(-9999,9999)
        charge.setValue(p.charge)
        form.addRow("Net molecular charge (e)",charge)
        if p.structure.name.lower().endswith(".pdb"):
            form.addRow(QLabel("Set the PDB charge for your chosen protonation state. The initial value is zero."))
        form.addRow(QLabel("Rotations use the laboratory X, Y and Z axes in order, then translation."))
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("OK")
        buttons.button(QDialogButtonBox.StandardButton.Cancel).setText("Cancel")
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        form.addRow(buttons)
        if not dlg.exec():
            return
        angles,shift = [tuple(spin.value() for spin in group) for group in fields]
        placement = Placement(p.structure,shift,angles,center.isChecked(),charge.value())
        label = (f"{p.structure.name}   ·   {len(p.structure.atoms):,} atoms   ·   "
                 f"XYZ {shift} Å   ·   rotation {angles}°   ·   charge {placement.charge:+d} e")
        if row is None:
            self.placements.append(placement)
            self.list.addItem(label)
            self.list.setCurrentRow(len(self.placements)-1)
        else:
            self.placements[row] = placement
            self.list.item(row).setText(label)
        self.changed.emit()
