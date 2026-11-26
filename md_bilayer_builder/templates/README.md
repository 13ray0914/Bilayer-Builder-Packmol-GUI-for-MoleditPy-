# Lipids of the Membrane recipe

One conformation of each of seven common phospholipids (CHARMM36 atom names and order),
a water molecule (TIP3P), and sodium and chloride ions: the built-in molecules of
`Membrane` and `MembraneProtein` (`src/Recipes/Membrane.jl`, `MembraneProtein.jl`). The lipids are the `<lipid>.pdb` files of the `packmol_lipids`
set (release `packmol_lipids-v1`) of [MolSimTestFiles](https://github.com/m3g/MolSimTestFiles),
whose README describes how they were obtained: the most rod-like frames of 2 ns of
restrained Langevin dynamics of each lipid in vacuum, centered at the origin, with the axis
from the midpoint of the chain ends to the head atom along +z. The same set has 20
conformations of each lipid (`<lipid>_conformations.pdb`).

| Lipid | Head atom | Chain-end atoms | Area per lipid (Å²) |
|:------|:---------:|:---------------:|:-------------------:|
| POPC  | 1 (N)     | 88 (C218), 131 (C316) | 64 |
| POPE  | 1 (N)     | 79 (C218), 122 (C316) | 57 |
| POPS  | 1 (N)     | 81 (C218), 124 (C316) | 63 |
| POPG  | 1 (C13)   | 81 (C218), 124 (C316) | 66 |
| DPPC  | 1 (N)     | 84 (C216), 127 (C316) | 63 |
| DOPC  | 1 (N)     | 88 (C218), 135 (C318) | 67 |
| DMPC  | 1 (N)     | 78 (C214), 115 (C314) | 61 |

The areas are those of the actual membranes, the defaults of `Membrane`.

## Additional PACKMOL-Memgen templates

| Family | Included species |
| --- | --- |
| Sterol | CHL1 (cholesterol) |
| Phosphatidylcholine | DLPC, DSPC, SDPC, DAPC |
| Phosphatidylethanolamine | DMPE, DPPE, DSPE, DOPE, DAPE |
| Phosphatidylserine | DMPS, DPPS, DOPS |
| Phosphatidylglycerol | DLPG, DMPG, DPPG, DSPG, DOPG |
| Phosphatidic acid | DPPA, POPA, DOPA |

The 21 templates in `memgen/` come from the PACKMOL-Memgen data distributed in `packmol_memgen_minimal` at commit `16ce513a677d808c438cb6d93b55fde6d7fdd945`. The PDB element fields were filled for packing and display. Coordinates and atom order are preserved. The accompanying `lipids.json` records one-based head and terminal atom indices, charge and reference area in Å². Areas use the experimental value from the parameter table when available, otherwise the force-field value. The polar oxygen is the head marker for cholesterol. The source license and attribution are included in `memgen/LICENSE` and `memgen/COPYRIGHT`. [Source templates and parameters](https://github.com/xiaoxuan-yu/packmol_memgen_minimal/tree/16ce513a677d808c438cb6d93b55fde6d7fdd945/packmol_memgen/data)
