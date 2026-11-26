# Bilayer Builder (Packmol-GUI for MoleditPy)

[![License](https://img.shields.io/github/license/13ray0914/Bilayer-Builder-Packmol-GUI-for-MoleditPy-?label=License&color=007ec6)](LICENSE)
[![Latest release](https://img.shields.io/github/v/release/13ray0914/Bilayer-Builder-Packmol-GUI-for-MoleditPy-?label=release&color=8957e5)](https://github.com/13ray0914/Bilayer-Builder-Packmol-GUI-for-MoleditPy-/releases/latest)

This independently developed MoleditPy plugin builds lipid bilayers with Packmol. Embed proteins and small molecules and inspect the packed structure in 3D.

## Installation

Download [Bilayer Builder (Packmol-GUI for MoleditPy).zip](https://github.com/13ray0914/Bilayer-Builder-Packmol-GUI-for-MoleditPy-/releases/latest/download/Bilayer.Builder.Packmol-GUI.for.MoleditPy.zip). Open **Plugin → Plugin Manager…** in MoleditPy and drag the ZIP into the manager. In the builder, select your Packmol executable under **Packmol Settings**.

See the [user guide](md_bilayer_builder/README.md) for leaflet composition, hydration, salt, embedded structures and 3D visualization. The [lipid catalog](md_bilayer_builder/templates/README.md) describes the included templates and their sources.

## License

The original Python code and documentation are available under the [MIT License](LICENSE). The bundled molecular templates and Windows runtime libraries have their own licenses. See the [third-party notices and corresponding source downloads](md_bilayer_builder/THIRD_PARTY.md).

## References

The bilayer setup draws on the membrane recipes in Packmol.jl. The PACKMOL-GUI paper informed the interface for configuring and previewing packing constraints.

1. *Packmol.jl Documentation*. [Stable documentation](https://m3g.github.io/Packmol.jl/stable/).
2. Huang, J., Wu, C., Yang, X., Yang, Z., Liu, S., and Yu, G. (2025). [PACKMOL-GUI: An All-In-One VMD Interface for Efficient Molecular Packing](https://pubs.acs.org/jcisd8/article/65/2/778/3686369/PACKMOL-GUI-An-All-In-One-VMD-Interface-for). *Journal of Chemical Information and Modeling*, 65(2), 778–784. [DOI 10.1021/acs.jcim.4c01639](https://doi.org/10.1021/acs.jcim.4c01639).
