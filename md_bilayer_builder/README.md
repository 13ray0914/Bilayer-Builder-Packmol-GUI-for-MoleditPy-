# Bilayer Builder (Packmol-GUI for MoleditPy)

Build a planar lipid bilayer in MoleditPy 5 with Packmol. Set the lipid composition of each leaflet independently, add TIP3P water and salt, and inspect the packed coordinates in 3D. The membrane lies in the XY plane, with its normal along Z.

## Installation

1. Open **Plugin → Plugin Manager…** in MoleditPy and drag the plugin ZIP into the manager. For a manual installation, copy the `md_bilayer_builder` folder into `~/.moleditpy/plugins/` and select **Reload Plugins**.
2. Open **Plugin → Bilayer Builder (Packmol-GUI for MoleditPy)**. Under **Packmol Settings**, select your Packmol executable and an existing output parent folder. Use a Packmol version that supports periodic boundaries. The plugin has been tested with MoleditPy 5.3.2 and Packmol 21.2.3 on Windows x64.
3. On Windows, the plugin looks first for `%LOCALAPPDATA%\Programs\Python\Python314\Scripts\packmol.exe`. Otherwise it looks for `packmol` on your PATH. You can browse to a different executable. The plugin remembers the full path and output location you select.

### Packmol 100.1.0

In **Packmol Settings**, clear **Start lipids and solvent on a grid**, then select **Build Bilayer**. This lets Packmol generate the initial positions. Clear the checkbox again after using **Restore Packing Defaults**.

If the **Run Log** shows `ERROR: restart_from placement.dat` (exit code 171), rerun the build with the checkbox cleared. The builder’s grid file stores molecule positions and rotations in the Fortran Packmol 21.2.3 format. Packmol 100.1.0 expects atomic coordinates (PDB, XYZ or mmCIF) for a restart. The `packall` and `movebadrandom` warnings can be ignored in 100.1.0. [Packmol migration guide](https://m3g.github.io/Packmol.jl/stable/migrating/)

## Build a bilayer

1. In **Lipid Composition**, expand **Sterols** or a phospholipid family (PA, PC, PE, PG or PS). PC opens by default. Each row shows the lipid type, charge, tail info, upper and lower leaflet counts, and reference area per lipid. **Tail Info** gives the carbon count and number of double bonds for sn-1 and sn-2, respectively. For example, POPC displays `16:0 / 18:1`. Cholesterol displays `N/A` because it has no sn-1/sn-2 acyl pair. Enter molecule counts for the upper (+Z) and lower (−Z) leaflets. Counts in collapsed groups are included in builds and in **Copy Upper Leaflet to Lower**.
2. To use a custom lipid, select **Add Lipid from PDB…** and choose a single-model PDB containing one lipid residue. Choose a family or use **Custom Lipids**, and enter tail information if known. Set the head atom index and the index of each tail end. Indices start at 1 and follow the order of atom records in the file. Enter the reference area per lipid in Å² and the molecule’s net charge in elementary charge units.
3. In **Cell, Water and Salt**, set the XY dimensions, water thickness on each side, salt type (NaCl or KCl), concentration and charge neutralization. The summary shows the planned molecule counts, net charge and added salt concentration after rounding to whole ion pairs.
4. To embed a protein or a synthetic molecule, open **Embedded Molecules → Add Structure…**. Select a single-model PDB (multiple chains and residues are supported) or a single-molecule 3D SDF/MOL file. Supply the protonation state and explicit atoms you want in the output. Enter the net charge for a PDB (default 0). SDF/MOL inputs start with their formal charge.
5. Choose the position and orientation of each structure. By default its geometric center moves to the origin, at the bilayer midplane. For a membrane-aligned PDB already in the desired coordinate frame, clear the centering box and leave all transformations at zero. Rotation around the laboratory X, Y and Z axes is applied in that order, followed by translation in Å. Use **Preview Constraints** to see the embedded atoms and the leaflet boundaries, and **Edit Placement…** to adjust the pose. Repeat **Add Structure…** for additional molecules or copies. Open **Packmol Input** to inspect the generated constraints.
6. Select **Build Bilayer** and follow the **Run Log**. Packmol holds the embedded structures fixed while placing the mobile molecules. The builder checks the final atom count and the identities and coordinates of the fixed atoms, then opens the 3D result. Drag to rotate, scroll to zoom, or select **Side View** or **Top View**. The **Embedded molecules** checkbox controls these atoms independently of the lipids, water and ions.

In **Packmol Settings**, enter an **Output folder name** such as `POPC-run-01`. Both **Save Input Files** and **Build Bilayer** create a new folder under the selected parent. If the name is already in use, the next export or build uses `POPC-run-01-2`, then `POPC-run-01-3`, and so on. Leave the name blank for an automatically generated `bilayer-…` folder. Select **Open Output Folder** to open the most recent export or build.

The **Recommended** column shows the builder’s starting values for bilayers. **Restore Packing Defaults** reapplies grid placement, a 2.00 Å heavy-atom tolerance, a 0.70 hydrogen radius scale, seed 12345, 200 maximum loops and 200 iterations per loop. The button applies to the six packing parameters.

The **Run Log** updates as Packmol writes output. During a long optimization step, the status line shows the elapsed run time and the time since the last output. A packing bar describes the current step. The completion message appears after the output atom count has been checked.

## Output files

| File | Contents |
| --- | --- |
| `bilayer.pdb` | Packed all-atom coordinates. The CRYST1 record gives the periodic cell. Chain IDs U and L identify the upper and lower leaflets. A and B identify the upper and lower solvent layers. |
| `packmol.inp` and component PDB files | Packing constraints and single-molecule templates. Run Packmol from this folder so the relative paths resolve. |
| `embedded_1.pdb`, `embedded_2.pdb`, … | Each imported structure at the requested pose. Protein chains and residue numbering are retained in the packed PDB. An embedded-atom count is recorded in the verified result for the 3D viewer. |
| `placement.dat` | Initial molecule positions and rotations when grid placement is enabled. Packmol optimizes both during packing. |
| `system.json` | Lipids, leaflet counts, reference areas, embedded structures and their poses and charges, settings, cell bounds, expected atom count and executable used. |
| `packmol.log` | Packmol output, errors and convergence information. |
| `bilayer.incomplete.pdb` | Intermediate coordinates, if Packmol wrote any before the run stopped or failed to converge. |

## Geometry and composition

- **Cell area.** For each leaflet, the builder sums the molecule count times the reference area of each lipid. It multiplies the larger leaflet area by the initial area scale factor (default 1.20) to size a square cell. Enter X and Y manually to use a rectangular cell. Both leaflets share the same area.
- **Head and tail constraints.** Heads face the water and tail ends lie near the midplane. The distance from the midplane to the head boundary is half the tail boundary distance plus the head-to-tail distance times the Z projection fraction. The head-to-tail distance is measured from the head atom to the mean position of the selected tail ends. On each side, the water layer starts 2 Å beyond the outermost head boundary. The cell is periodic in X, Y and Z.
- **Initial placement.** By default, lipids and solvent start on a grid. The random seed controls lipid assignment to grid points and in-plane rotation. Clear **Start lipids and solvent on a grid** to let Packmol generate the starting positions.
- **Embedded structures.** With automatic sizing, the builder adds the estimated cross-sectional area of the embedded atoms to each leaflet's target area. A 2 Å voxel grid estimates the footprint and displaced solvent. The summary uses the remaining water-layer volume to calculate hydration and salt. Fixed molecules must fit inside the cell and be separated from each other by the packing tolerance.
- **Salt and counterions.** The number of salt pairs is the requested molarity times the accessible water volume in liters times the Avogadro constant, rounded to the nearest whole pair. The displayed added salt concentration uses this integer count. Charge neutralization adds counterions separately to balance the combined charge of the lipids and embedded structures.
- **Water count.** Under **Cell, Water and Salt → Water count method**, choose **Estimate from water volume** to start at 1 g/cm³ and replace one water with each ion. **Water molecules per lipid** multiplies the ratio by the total lipid count in both leaflets (including sterols) and rounds to the nearest whole water, with halves rounded up. For 64 upper and 64 lower lipids, a ratio of 30 gives 3,840 waters. **Total water molecules** accepts an exact integer for the whole cell. Both explicit methods preserve the requested number when adding ions and distribute waters in proportion to the accessible volume on each side. Water thickness sets the available packing space and affects the salt calculation in all modes. Increase it if the specified waters need more space. For a dry bilayer, set both water thickness and salt concentration to zero and disable neutralization if the system carries a net charge.
- **Packing radii.** The heavy-atom radius is half the Packmol tolerance. The default hydrogen radius is 0.70 times the heavy-atom radius. Set the hydrogen scale to 1.00 to apply the same distance criterion to all atoms. If a mixture fails to converge, inspect the run log and constraint preview before adjusting the XY area or head-to-tail projection.
- **Interpreting the structure.** The viewer infers backbone lines from heavy-atom distances within each lipid residue. PDB output uses the CHARMM convention for four-character residue names. Before molecular dynamics, match the coordinates to the lipid, water and ion topologies for your force field, minimize the energy and equilibrate the system.

## Sources and runtime

The original Python code and documentation are distributed under the [MIT License](LICENSE). See the [third-party notices](THIRD_PARTY.md) for the licenses and corresponding source of the bundled materials.

The library includes 28 lipid templates, covering cholesterol (CHL1) and phosphatidylcholine, phosphatidylethanolamine, phosphatidylserine, phosphatidylglycerol and phosphatidic acid species. The seven original lipid templates come from Packmol.jl and the other 21 from PACKMOL-Memgen data. See [the template catalog](templates/README.md) for names, provenance and licenses. [Packmol.jl](https://github.com/m3g/Packmol.jl/tree/f2d7af921be07b7ac9ec14279a47e842247313ca/src/Recipes/lipids) · [PACKMOL-Memgen](https://doi.org/10.1021/acs.jcim.9b00269)

The leaflet controls draw on CHARMM-GUI Membrane Builder. The constraint preview, input inspection and run log draw on PACKMOL-GUI. [CHARMM-GUI Membrane Builder](https://www.charmm-gui.org/?doc=demo&id=membrane_builder&lesson=3) · [PACKMOL-GUI paper](https://doi.org/10.1021/acs.jcim.4c01639)

Packmol runs as an external program. For the Windows x64 pip distribution of Packmol 21.2.3, the plugin includes MSYS2 UCRT64 runtime DLLs from gcc-libgfortran and gcc-libs 15.2.0-9 and libwinpthread 13.0.0.r566.g9b3dd0125-1. The runtime search path applies to the launched Packmol process. License files are in `bin/windows-x64/runtime-licenses/`. [Packmol](https://github.com/m3g/packmol) · [MSYS2 packages](https://repo.msys2.org/mingw/ucrt64/)
