# Third-party notices

## Molecular templates

- The seven Packmol.jl lipids and its water, sodium and chloride templates are distributed with the upstream [MIT notice](templates/LICENSE). The [template catalog](templates/README.md) lists their origins and atom indices.
- The 21 additional PDB templates and the parameters derived from PACKMOL-Memgen use GPL-2.0-only. The corresponding editable PDB and JSON files, [license](templates/memgen/LICENSE) and [attribution](templates/memgen/COPYRIGHT) are included in the plugin. Element fields were filled in for packing and display. The source data and parameter table are available at the [upstream revision](https://github.com/xiaoxuan-yu/packmol_memgen_minimal/tree/16ce513a677d808c438cb6d93b55fde6d7fdd945/packmol_memgen/data).

## Windows runtime libraries

The Windows x64 plugin includes four unmodified MSYS2 UCRT64 DLLs for the external Packmol process.

| Libraries | Distribution terms and notices | Corresponding source |
| --- | --- | --- |
| `libgcc_s_seh-1.dll`, `libgfortran-5.dll` | GPL-3.0-or-later with the GCC Runtime Library Exception 3.1. [GPL](bin/windows-x64/runtime-licenses/gcc-libs/COPYING3) · [Runtime exception](bin/windows-x64/runtime-licenses/gcc-libgfortran/COPYING.RUNTIME) | [MSYS2 GCC 15.2.0-9 source and packaging patches](https://github.com/13ray0914/Bilayer-Builder-Packmol-GUI-for-MoleditPy-/releases/download/v1.0.1/mingw-w64-gcc-15.2.0-9.src.tar.zst) |
| `libquadmath-0.dll` | LGPL-2.1-or-later. [License](bin/windows-x64/runtime-licenses/gcc-libs/COPYING.LIB) · [MSYS2 notice](bin/windows-x64/runtime-licenses/gcc-libs/README) | Included in the same GCC source archive. [Upstream mirror](https://repo.msys2.org/mingw/sources/mingw-w64-gcc-15.2.0-9.src.tar.zst) |
| `libwinpthread-1.dll` | The mingw-w64 project and Lockless Inc. notices apply to their respective portions. [License and copyright notices](bin/windows-x64/runtime-licenses/libwinpthread/COPYING) | [MSYS2 winpthreads 13.0.0.r566.g9b3dd0125-1 source](https://repo.msys2.org/mingw/sources/mingw-w64-winpthreads-13.0.0.r566.g9b3dd0125-1.src.tar.zst) |
