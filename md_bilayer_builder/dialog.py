import codecs
import io
import json
import os
import platform
import shutil
import subprocess
from contextlib import ExitStack
from pathlib import Path

from PySide6.QtCore import QElapsedTimer, QProcess, QProcessEnvironment, QSettings, QSignalBlocker, QTimer, Qt, QUrl
from PySide6.QtGui import QDesktopServices, QFontDatabase, QTextCursor
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDialog, QDialogButtonBox,
    QDoubleSpinBox, QFileDialog, QFormLayout, QGridLayout, QHBoxLayout, QHeaderView, QLabel, QLineEdit,
    QMessageBox, QPlainTextEdit, QPushButton, QSpinBox, QTabWidget, QTreeWidget,
    QTreeWidgetItem, QVBoxLayout, QWidget)

from . import PLUGIN_NAME
from .model import LIPID_GROUPS, Lipid, Plan, Settings, library, molecule_records
from .embedded_panel import EmbeddedPanel
from .inserts import verify_fixed_output


def number(value, lo, hi, decimals=0):
    if decimals:
        field = QDoubleSpinBox()
        field.setDecimals(decimals)
        field.setSingleStep(0.01 if decimals >= 2 else 0.5)
    else:
        field = QSpinBox()
    field.setRange(lo, hi)
    field.setValue(value)
    field.setKeyboardTracking(False)
    return field


def button(text, callback):
    btn = QPushButton(text)
    btn.clicked.connect(callback)
    return btn


def editor():
    field = QPlainTextEdit()
    field.setReadOnly(True)
    field.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
    field.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
    return field


def default_executable():
    if os.name == "nt":
        local = os.environ.get("LOCALAPPDATA", "")
        candidate = Path(local) / "Programs/Python/Python314/Scripts/packmol.exe"
        if candidate.is_file():
            return str(candidate)
    return shutil.which("packmol") or ""


class Builder(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(PLUGIN_NAME)
        self.resize(810, 655)
        # Keep the settings namespace used by earlier releases.
        self.prefs = QSettings("MD Bilayer Builder", "Packmol")
        self.lipid_rows = []
        self.proc = None
        self.folder = self.result = None
        self.viewers = []
        root = QVBoxLayout(self)
        self.tabs = QTabWidget()
        root.addWidget(self.tabs, 1)
        lipids = QWidget()
        self.tabs.addTab(lipids, "Lipid Composition")
        layout = QVBoxLayout(lipids)
        layout.addWidget(QLabel("Enter the number of molecules of each lipid in each leaflet. The membrane normal is along the Z axis."))
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Lipid Type", "Charge\n(e)", "Tail Info\n(sn-1 / sn-2)",
                                   "Upper Leaflet\n(+Z)", "Lower Leaflet\n(−Z)", "Area per Lipid\n(Å²)"])
        self.tree.headerItem().setToolTip(2, "Carbon atoms and double bonds in each acyl chain, in sn-1 / sn-2 order.")
        header = self.tree.header()
        header.setStretchLastSection(False)
        header.setSectionsMovable(False)
        header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        for column in (3, 4):
            header.setSectionResizeMode(column, QHeaderView.ResizeMode.Stretch)
        self.tree.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        layout.addWidget(self.tree)
        self.groups = {}
        for key, title in LIPID_GROUPS.items():
            group = QTreeWidgetItem(self.tree, [title])
            group.setFirstColumnSpanned(True)
            font = group.font(0)
            font.setBold(True)
            group.setFont(0, font)
            group.setHidden(True)
            self.groups[key] = group
        for lipid in library():
            self.add_row(lipid, 16 if lipid.name == "POPC" else 0)
        self.groups["PC"].setExpanded(True)
        actions = QHBoxLayout()
        layout.addLayout(actions)
        actions.addWidget(button("Add Lipid from PDB…", self.import_lipid))
        actions.addWidget(button("Copy Upper Leaflet to Lower", self.copy_upper))
        actions.addStretch()

        box = QWidget()
        self.tabs.addTab(box, "Cell, Water and Salt")
        form = QFormLayout(box)
        self.auto_xy = QCheckBox("Calculate from the leaflet compositions")
        self.auto_xy.setChecked(True)
        form.addRow("Cell dimensions in XY", self.auto_xy)
        self.scale = number(1.20, .50, 5, 2)
        form.addRow("Initial area scale factor", self.scale)
        self.x, self.y = number(50., 5., 2000., 1), number(50., 5., 2000., 1)
        form.addRow("Manual X length (Å)", self.x)
        form.addRow("Manual Y length (Å)", self.y)
        self.water = number(18., 0., 300., 1)
        form.addRow("Water thickness on each side (Å)", self.water)
        self.water_mode = QComboBox()
        self.water_mode.addItem("Estimate from water volume", "volume")
        self.water_mode.addItem("Water molecules per lipid", "per_lipid")
        self.water_mode.addItem("Total water molecules", "total")
        form.addRow("Water count method", self.water_mode)
        self.waters_per_lipid = number(Settings.waters_per_lipid, 0, 10000, 2)
        self.total_waters = number(Settings.total_waters, 0, 2_000_000_000)
        self.waters_per_lipid.setToolTip("Multiply by the total number of lipids across both leaflets.")
        self.water.setToolTip("Sets the packing space and the volume used for salt concentration in all three methods.")
        form.addRow("Water molecules per lipid", self.waters_per_lipid)
        form.addRow("Total water molecules", self.total_waters)
        self.salt_type = QComboBox()
        self.salt_type.addItems(["NaCl", "KCl"])
        form.addRow("Salt", self.salt_type)
        self.salt = number(.15, 0., 5., 3)
        form.addRow("Added salt in the water layers (mol/L)", self.salt)
        self.neutral = QCheckBox("Add counterions to neutralize the system")
        self.neutral.setChecked(True)
        form.addRow("Charge neutralization", self.neutral)
        self.tilt = number(.88, .25, 1., 2)
        self.tail = number(4., 1., 10., 1)
        form.addRow("Head-to-tail projection on Z (fraction)", self.tilt)
        form.addRow("Tail boundary from the midplane (Å)", self.tail)

        self.embedded = EmbeddedPanel()
        self.tabs.addTab(self.embedded, "Embedded Molecules")
        engine = QWidget()
        self.tabs.addTab(engine, "Packmol Settings")
        form = QFormLayout(engine)
        self.exe = QLineEdit(self.prefs.value("exe", default_executable(), type=str))
        self.output = QLineEdit(self.prefs.value("output", str(Path.home()), type=str))
        for label, field, browse in [("Packmol executable", self.exe, self.choose_exe),
                                     ("Output parent folder", self.output, self.choose_output)]:
            row = QHBoxLayout()
            row.addWidget(field)
            row.addWidget(button("Browse…", browse))
            form.addRow(label, row)
        self.folder_name = QLineEdit()
        self.folder_name.setPlaceholderText("Leave blank for an automatic name")
        form.addRow("Output folder name", self.folder_name)
        form.addRow("", QLabel("Existing names get a numbered suffix, starting with -2."))

        defaults = Settings()
        self.tolerance = number(defaults.tolerance, .5, 5., 2)
        self.hscale = number(defaults.hydrogen_scale, .50, 1., 2)
        self.seed = number(defaults.seed, 1, 2_000_000_000)
        self.loops = number(defaults.loops, 1, 5000)
        self.maxit = number(defaults.maxit, 1, 5000)
        self.lattice = QCheckBox("Start lipids and solvent on a grid")
        self.lattice.setChecked(defaults.lattice)
        self.packing_defaults = (
            ("Initial placement", self.lattice, defaults.lattice),
            ("Heavy-atom tolerance (Å)", self.tolerance, defaults.tolerance),
            ("Hydrogen packing radius scale", self.hscale, defaults.hydrogen_scale),
            ("Random seed", self.seed, defaults.seed),
            ("Maximum loops (nloop)", self.loops, defaults.loops),
            ("Iterations per loop (maxit)", self.maxit, defaults.maxit),
        )
        grid = QGridLayout()
        grid.setColumnStretch(1, 1)
        for col, title in enumerate(("Packing parameter", "Value", "Recommended")):
            grid.addWidget(QLabel(title), 0, col)
        for row, (label, field, value) in enumerate(self.packing_defaults, 1):
            caption = QLabel(label)
            caption.setBuddy(field)
            grid.addWidget(caption, row, 0)
            grid.addWidget(field, row, 1)
            if field is self.lattice:
                recommended = "On" if value else "Off"
            else:
                recommended = field.textFromValue(value)
            grid.addWidget(QLabel(recommended), row, 2)
        form.addRow(grid)
        form.addRow(button("Restore Packing Defaults", self.restore_packing_defaults))
        self.input = editor()
        self.tabs.addTab(self.input, "Packmol Input")
        self.log = editor()
        self.log.document().setMaximumBlockCount(600)
        self.tabs.addTab(self.log, "Run Log")
        self.summary = QLabel()
        self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        root.addWidget(self.summary)
        bar = QHBoxLayout()
        root.addLayout(bar)
        self.preview_btn = button("Preview Constraints", self.preview)
        self.export_btn = button("Save Input Files", self.export_only)
        self.run_btn = button("Build Bilayer", self.start)
        self.stop_btn = button("Stop", self.stop)
        self.open_btn = button("View Result in 3D…", self.open_result)
        self.dir_btn = button("Open Output Folder", self.open_folder)
        for b in (self.preview_btn, self.export_btn, self.run_btn, self.stop_btn, self.open_btn, self.dir_btn):
            bar.addWidget(b)
        self.status = QLabel("Ready")
        self.status.setWordWrap(True)
        self.status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        root.addWidget(self.status)
        self.timer = QTimer(self)
        self.timer.setInterval(500)
        self.timer.timeout.connect(self.read_log)
        for field in [self.scale, self.x, self.y, self.water, self.salt, self.tail, self.tilt,
                      self.tolerance, self.hscale, self.seed, self.loops, self.maxit,
                      self.waters_per_lipid, self.total_waters]:
            field.valueChanged.connect(self.refresh)
        self.salt_type.currentIndexChanged.connect(self.refresh)
        self.water_mode.currentIndexChanged.connect(self.refresh)
        self.neutral.toggled.connect(self.refresh)
        self.auto_xy.toggled.connect(self.refresh)
        self.lattice.toggled.connect(self.refresh)
        self.embedded.changed.connect(self.refresh)
        self.refresh()
        from PySide6.QtWidgets import QApplication
        QApplication.instance().aboutToQuit.connect(self.stop)

    def add_row(self, lipid, count):
        group = self.groups[lipid.group]
        group.setHidden(False)
        label = "CHL1 (cholesterol)" if lipid.name == "CHL1" and lipid.group == "Sterols" else lipid.name
        item = QTreeWidgetItem(group, [label, str(lipid.charge), lipid.tail_info, "", "", f"{lipid.area:g}"])
        for column in (1, 2, 5):
            item.setTextAlignment(column, Qt.AlignmentFlag.AlignCenter)
        upper, lower = number(count, 0, 9999), number(count, 0, 9999)
        for column, field, side in ((3, upper, "upper"), (4, lower, "lower")):
            field.setAccessibleName(f"{lipid.name} {side} leaflet count")
            self.tree.setItemWidget(item, column, field)
            field.valueChanged.connect(self.refresh)
        self.lipid_rows.append((lipid, upper, lower))
        return item

    def copy_upper(self):
        with ExitStack() as signals:
            for _, upper, lower in self.lipid_rows:
                signals.enter_context(QSignalBlocker(lower))
                lower.setValue(upper.value())
        self.refresh()

    def import_lipid(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open a single-lipid PDB", "", "PDB (*.pdb);;All files (*)")
        if not path:
            return
        try:
            atoms = molecule_records(path)
            dlg = QDialog(self)
            dlg.setWindowTitle(f"Lipid Template  ·  {Path(path).name}")
            form = QFormLayout(dlg)
            labels = editor()
            labels.setPlainText("\n".join(f"{i:4}  {s[12:16].strip()}  {s[17:21].strip()}" for i,s in enumerate(atoms,1)))
            form.addRow("Atom index / name / residue", labels)
            group = QComboBox()
            for key, title in LIPID_GROUPS.items():
                group.addItem(title, key)
            group.setCurrentIndex(group.findData("Custom"))
            tail_info = QLineEdit()
            tail_info.setPlaceholderText("Optional, e.g. 16:0 / 18:1")
            form.addRow("Lipid group", group)
            form.addRow("Tail info (sn-1 / sn-2)", tail_info)
            head, tails = number(1, 1, len(atoms)), QLineEdit()
            tails.setPlaceholderText("e.g., 88, 131")
            area, charge = number(64., 1., 1000., 1), number(0,-20,20)
            for name, field in [("Head atom index",head),("Tail-end atom indices",tails),("Area per lipid (Å²)",area),("Net charge (e)",charge)]:
                form.addRow(name, field)
            buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
            buttons.button(QDialogButtonBox.StandardButton.Ok).setText("OK")
            buttons.button(QDialogButtonBox.StandardButton.Cancel).setText("Cancel")
            buttons.accepted.connect(dlg.accept)
            buttons.rejected.connect(dlg.reject)
            form.addRow(buttons)
            if dlg.exec():
                ends = tuple(int(x) for x in tails.text().replace(",", " ").split())
                lipid = Lipid(Path(path).stem, atoms, head.value(), ends, area.value(), charge.value(),
                              group.currentData(), tail_info.text().strip() or "Unspecified")
                item = self.add_row(lipid, 0)
                item.parent().setExpanded(True)
                self.tree.setCurrentItem(item)
                self.tree.scrollToItem(item)
                self.refresh()
        except (OSError, ValueError) as exc:
            self.error(exc)

    def plan(self):
        rows = [(m, upper.value(), lower.value()) for m, upper, lower in self.lipid_rows]
        settings = Settings(water=self.water.value(), salt=self.salt.value(), cation=["Na","K"][self.salt_type.currentIndex()],
                            neutralize=self.neutral.isChecked(), xy=() if self.auto_xy.isChecked() else (self.x.value(), self.y.value()),
                            area_scale=self.scale.value(), extension=self.tilt.value(), tail_band=self.tail.value(),
                            tolerance=self.tolerance.value(), hydrogen_scale=self.hscale.value(),
                            seed=self.seed.value(), loops=self.loops.value(), maxit=self.maxit.value(),
                            lattice=self.lattice.isChecked(), water_mode=self.water_mode.currentData(),
                            waters_per_lipid=self.waters_per_lipid.value(), total_waters=self.total_waters.value())
        return Plan(rows, settings, self.embedded.placements)

    def refresh(self):
        auto = self.auto_xy.isChecked()
        self.x.setEnabled(not auto)
        self.y.setEnabled(not auto)
        self.scale.setEnabled(auto)
        self.waters_per_lipid.setEnabled(self.water_mode.currentData() == "per_lipid")
        self.total_waters.setEnabled(self.water_mode.currentData() == "total")
        try:
            plan = self.plan()
            self.summary.setText(plan.summary)
            self.input.setPlainText(plan.deck())
            valid = True
        except ValueError as exc:
            self.summary.setText(str(exc))
            self.input.clear()
            valid = False
        for btn in (self.preview_btn, self.run_btn, self.export_btn):
            btn.setEnabled(valid and self.proc is None)
        self.stop_btn.setEnabled(self.proc is not None)
        self.dir_btn.setEnabled(self.folder is not None)
        for i in range(4):
            self.tabs.setTabEnabled(i, self.proc is None)

    def choose_exe(self):
        path, _ = QFileDialog.getOpenFileName(self, "Packmol executable", self.exe.text())
        if path:
            self.exe.setText(path)

    def choose_output(self):
        path = QFileDialog.getExistingDirectory(self, "Select Output Folder", self.output.text())
        if path:
            self.output.setText(path)

    def restore_packing_defaults(self):
        with ExitStack() as signals:
            for _, field, value in self.packing_defaults:
                signals.enter_context(QSignalBlocker(field))
                if field is self.lattice:
                    field.setChecked(value)
                else:
                    field.setValue(value)
        self.refresh()

    def error(self, exc):
        self.status.setText(str(exc))
        QMessageBox.warning(self, PLUGIN_NAME, str(exc))

    def export(self, plan, executable=None):
        parent = Path(self.output.text().strip()).expanduser()
        if not self.output.text().strip() or not parent.is_dir():
            raise ValueError("Choose an existing output folder in Packmol Settings.")
        self.folder = plan.export(parent.resolve(), executable, folder_name=self.folder_name.text())
        self.prefs.setValue("output", str(parent))
        self.refresh()

    def export_only(self):
        try:
            self.export(self.plan())
            self.status.setText(f"Input files saved  ·  {self.folder}")
        except (ValueError, OSError) as exc:
            self.error(exc)

    def start(self):
        if self.proc:
            return
        try:
            exe = Path(self.exe.text().strip().strip('"')).expanduser()
            if not exe.is_absolute() or not exe.is_file():
                raise ValueError("Enter the full path to the Packmol executable in Packmol Settings.")
            self.active_plan = self.plan()
            self.export(self.active_plan, exe)
        except (ValueError, OSError) as exc:
            self.error(exc)
            return
        self.prefs.setValue("exe", str(exe))
        self.log.clear()
        self.log.setPlainText(f"{exe}\n{self.folder}\n")
        self.result = None
        self.offset = 0
        self.log_decoder = io.IncrementalNewlineDecoder(
            codecs.getincrementaldecoder("utf-8")(errors="replace"), translate=True)
        self.elapsed = QElapsedTimer()
        self.elapsed.start()
        self.last_output_ms = 0
        self.cancelled = False
        proc = QProcess(self)
        self.proc = proc
        env = QProcessEnvironment.systemEnvironment()
        # Flush Fortran console output while keeping bulk PDB file writes buffered.
        env.insert("GFORTRAN_UNBUFFERED_PRECONNECTED", "1")
        runtime = Path(__file__).parent / "bin/windows-x64"
        # The pip launcher starts the actual Fortran binary as a child.
        # The Windows pip wheel omits its UCRT Fortran DLLs. Older MinGW DLLs
        # elsewhere on PATH can share a name but use an incompatible ABI.
        wheel_binary = exe.parent.parent / "Lib/site-packages/packmol/binaries/packmol.exe"
        if (os.name == "nt" and platform.machine().lower() in ("amd64","x86_64")
                and exe.parent.name.lower() == "scripts" and wheel_binary.is_file()):
            env.insert("PATH", str(runtime) + os.pathsep + env.value("PATH"))
        proc.setProcessEnvironment(env)
        proc.setWorkingDirectory(str(self.folder))
        proc.setStandardInputFile(str(self.folder / "packmol.inp"))
        proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        proc.setStandardOutputFile(str(self.folder / "packmol.log"))
        proc.finished.connect(self.packing_finished)
        proc.errorOccurred.connect(self.process_error)
        self.refresh()
        self.tabs.setCurrentWidget(self.log)
        self.status.setText(f"Running Packmol  ·  {self.folder}")
        self.timer.start()
        proc.start(str(exe), [])

    def read_log(self, final=False):
        path = self.folder / "packmol.log"
        data = b""
        if path.exists():
            with path.open("rb") as stream:
                stream.seek(self.offset)
                data = stream.read()
                self.offset = stream.tell()
        now = self.elapsed.elapsed()
        if data:
            self.last_output_ms = now
        chunk = self.log_decoder.decode(data, final=final)
        if chunk:
            scrollbar = self.log.verticalScrollBar()
            follow = scrollbar.value() == scrollbar.maximum()
            cursor = QTextCursor(self.log.document())
            cursor.movePosition(QTextCursor.MoveOperation.End)
            cursor.insertText(chunk)
            if follow:
                scrollbar.setValue(scrollbar.maximum())
        if not final:
            minutes, seconds = divmod(now // 1000, 60)
            quiet = (now - self.last_output_ms) // 1000
            activity = f"Last output {quiet}s ago" if self.offset else "Waiting for first output"
            self.status.setText(f"Running Packmol  ·  {minutes}:{seconds:02d} elapsed  ·  {activity}  ·  {self.folder}")

    def release_process(self):
        self.timer.stop()
        self.read_log(final=True)
        proc, self.proc = self.proc, None
        proc.deleteLater()
        self.refresh()

    def process_error(self, error):
        if self.proc and error == QProcess.ProcessError.FailedToStart:
            msg = self.proc.errorString()
            self.release_process()
            self.log.appendPlainText(f"Could not start Packmol. {msg}")
            self.error(f"Could not start Packmol. {msg}")

    def packing_finished(self, code, exit_status):
        if not self.proc:
            return
        self.release_process()
        path = self.folder / "bilayer.pdb"
        if self.cancelled:
            if path.is_file():
                path.rename(path.with_name("bilayer.incomplete.pdb"))
            self.status.setText(f"Run stopped  ·  {self.folder}")
            self.log.appendPlainText("Run stopped.")
            return
        try:
            with (self.folder / "packmol.log").open("rb") as stream:
                success = any(line.strip().lower() == b"success!" for line in stream)
            if exit_status != QProcess.ExitStatus.NormalExit or code or not success:
                raise ValueError(f"Packmol did not finish successfully (exit code {code}). See the Run Log for details.")
            with path.open("rb") as stream:
                total = sum(line.startswith((b"ATOM  ",b"HETATM")) for line in stream)
            if total != self.active_plan.atom_count:
                raise ValueError(f"Expected {self.active_plan.atom_count:,} atoms in the output PDB, but found {total:,}.")
            verify_fixed_output(path, self.active_plan.placements)
            self.result = path
            self.status.setText(f"Bilayer complete  ·  {total:,} atoms  ·  {path}")
            self.log.appendPlainText(self.status.text())
            self.show_view(path=path, bounds=self.active_plan.bounds)
        except (OSError, ValueError) as exc:
            if path.is_file():
                path.rename(path.with_name("bilayer.incomplete.pdb"))
            self.log.appendPlainText(str(exc))
            self.error(exc)

    def stop(self):
        if self.proc is None:
            return
        self.cancelled = True
        pid = self.proc.processId()
        if os.name == "nt" and pid:
            subprocess.run([str(Path(os.environ["SystemRoot"]) / "System32/taskkill.exe"), "/PID", str(pid), "/T", "/F"],
                           creationflags=subprocess.CREATE_NO_WINDOW, capture_output=True, check=False)
        if self.proc:
            self.proc.kill()

    def preview(self):
        self.show_view(plan=self.plan())

    def open_result(self):
        path = self.result
        if not path:
            chosen, _ = QFileDialog.getOpenFileName(self, "Open a packed bilayer", str(self.folder or ""), "PDB (*.pdb)")
            if not chosen:
                return
            path = Path(chosen)
        manifest = path.with_name("system.json")
        try:
            bounds = json.loads(manifest.read_text(encoding="utf-8")).get("bounds") if manifest.is_file() else None
            self.show_view(path=path, bounds=bounds)
        except (ValueError, OSError) as exc:
            self.error(exc)

    def show_view(self, **kwargs):
        try:
            from .viewer import Viewer
            win = Viewer(self, **kwargs)
            self.viewers.append(win)
            win.destroyed.connect(lambda *_: self.viewers.remove(win) if win in self.viewers else None)
            win.show()
        except Exception as exc:
            self.error(f"Could not open the 3D viewer. {exc}")

    def open_folder(self):
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.folder)))

    def reject(self):
        if self.proc is None:
            super().reject()
        else:
            self.status.setText("Packmol is running. Select Stop before closing this window.")
