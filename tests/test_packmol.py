"""Exercise the builder's normal Packmol launch and result handling."""

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QApplication

from md_bilayer_builder.dialog import Builder, default_executable, resolve_executable
from md_bilayer_builder.viewer import read_atoms


class PackmolTest(unittest.TestCase):
    def test_build_in_active_python_environment(self):
        app = QApplication.instance() or QApplication([])
        with tempfile.TemporaryDirectory(prefix="bilayer test ") as folder:
            QSettings.setDefaultFormat(QSettings.Format.IniFormat)
            for scope in (QSettings.Scope.UserScope, QSettings.Scope.SystemScope):
                QSettings.setPath(QSettings.Format.IniFormat, scope, folder)
            # The launcher installed with the active Python is also discoverable
            # when the GUI inherits a PATH without the Python scripts directory.
            with patch.dict(os.environ, {"PATH": ""}):
                executable = default_executable()
            self.assertTrue(executable, "Install Packmol 21.2.3 in this Python environment.")
            self.assertEqual(resolve_executable(f'"{executable}"'), Path(executable))
            builder = Builder()
            builder.output.setText(folder)
            builder.exe.setText(executable)
            builder.auto_xy.setChecked(False)
            builder.x.setValue(28)
            builder.y.setValue(32)
            builder.water.setValue(12)
            builder.salt.setValue(0)
            builder.water_mode.setCurrentIndex(builder.water_mode.findData("total"))
            builder.total_waters.setValue(12)
            for lipid, upper, lower in builder.lipid_rows:
                upper.setValue(2 if lipid.name == "POPC" else 0)
                lower.setValue(2 if lipid.name == "DMPC" else 0)
            # Keep this test independent of the runner's graphics hardware.
            with patch.object(builder, "error") as error, patch.object(builder, "show_view") as view:
                builder.start()
                process = builder.proc
                if process:
                    finished = process.waitForFinished(45_000)
                    if not finished and builder.proc:
                        builder.stop()
                        process.waitForFinished(5_000)
                    self.assertTrue(finished, builder.log.toPlainText())
                self.assertFalse(error.called, builder.log.toPlainText())
                self.assertIsNotNone(builder.result, builder.log.toPlainText())
                view.assert_called_once_with(path=builder.result, bounds=builder.active_plan.cell_bounds)
                _, bounds, _ = read_atoms(builder.result)
                self.assertEqual(bounds[::2], (0., 0., 0.))
                for actual, expected in zip(bounds, builder.active_plan.cell_bounds):
                    self.assertAlmostEqual(actual, expected, delta=.006)
            builder.close()
            app.processEvents()
