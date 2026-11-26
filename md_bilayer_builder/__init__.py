PLUGIN_NAME = "Bilayer Builder (Packmol-GUI for MoleditPy)"
PLUGIN_VERSION = "1.0.2"
PLUGIN_AUTHOR = "13ray0914"
PLUGIN_LICENSE = "MIT"
PLUGIN_DESCRIPTION = "Build lipid bilayers with Packmol, embed proteins and small molecules, and inspect the packed structure in 3D."
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=5.0.0, <6.0.0"
PLUGIN_SUPPORTED_PYTHON_VERSION = ">=3.9, <3.15"
PLUGIN_TAGS = ["Calculation", "Visualization"]
PLUGIN_DEPENDENCIES = ["PySide6", "numpy", "pyvista", "pyvistaqt"]
PLUGIN_OPTIONAL_DEPENDENCIES = ["rdkit"]


def initialize(context):
    def open_builder():
        dialog = context.get_window("builder")
        if dialog is None:
            from .dialog import Builder
            dialog = Builder(context.get_main_window())
            context.register_window("builder", dialog)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    context.add_menu_action(f"Plugin/{PLUGIN_NAME}", open_builder)
