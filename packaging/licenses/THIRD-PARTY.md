# PrimerPicker: licences

PrimerPicker is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version (`LICENSE`). It is distributed WITHOUT
ANY WARRANTY. Its source code is at
<https://github.com/frederic-mahe/primerpicker>, with a source archive
attached to each release.

The packages of PrimerPicker also contain software and data written by
others, under their own licences, listed below.

## REBASE

The restriction enzymes come from REBASE, The Restriction Enzyme
Database (<http://rebase.neb.com>), version 610, which may be
distributed with software provided users are not charged for the
REBASE data and the following reference is cited:

> Roberts, R.J., Vincze, T., Posfai, J., Macelis, D. (2010)
> REBASE--a database for DNA restriction and modification: enzymes,
> genes and genomes. Nucl. Acids Res. 38: D234-D236.

## Qt (the graphical program)

The graphical program uses the Qt toolkit (Copyright The Qt Company Ltd
and other contributors), under the GNU Lesser General Public License,
version 3 (`LGPL-3.0.txt`, which adds permissions to the GNU General
Public License, version 3, `LICENSE`). Qt is not modified, and it is
linked dynamically: its libraries are separate files of the package
(DLL files next to `primerpicker-gui.exe` on Windows, `usr/lib` inside
the AppImage on Linux, `Contents/Frameworks` inside `PrimerPicker.app`
on macOS), which can be replaced with other builds of the same Qt
version. Help -> About PrimerPicker shows the Qt version.

The source code of every Qt version is at
<https://download.qt.io/archive/qt/>, in the folder of that version
(for Qt 6.12.0: `6.12/6.12.0/single/`). Qt contains code from other
projects (for instance PCRE2, HarfBuzz, libpng, libjpeg, the Unicode
data), each under its own licence; the Linux and macOS packages
include the software bill of materials of Qt (`qt-sbom/`), which
lists them with their licences, and Qt's documentation lists them too:
<https://doc.qt.io/qt-6/licenses-used-in-qt.html>.

## Libraries of each package

- **Windows** (graphical program): Qt and the other libraries come
  from MSYS2 (<https://www.msys2.org>, UCRT64 environment). Each
  MSYS2 package with a library in the folder has its licence files in
  `licenses/<package>/`, and `licenses/packages.txt` gives the exact
  package versions; their sources are at
  <https://github.com/msys2/MINGW-packages> and
  <https://repo.msys2.org/mingw/sources/>.
- **Linux** (AppImage): Qt (official binaries) and ICU 73
  (`ICU-73.txt`, Unicode licence) come from the Qt Company; the other
  libraries come from Ubuntu 24.04, each with its copyright file in
  `usr/share/doc/<package>/copyright`, and their sources at
  <https://launchpad.net/ubuntu/noble>. The AppImage runtime (the
  start of the AppImage file, <https://github.com/AppImage/type2-runtime>,
  release 20251108) is under the MIT licence and contains libfuse
  (GNU Lesser General Public License, version 2.1, `LGPL-2.1.txt`),
  squashfuse (BSD 2-clause) and zstd (BSD 3-clause). The files of the
  AppImage are listed by `primerpicker-gui-*.AppImage --appimage-extract`.
- **macOS** (disk image): Qt (official binaries) only; the other
  libraries are those of macOS.
- The command-line programs contain the GCC runtime libraries (GNU
  General Public License, version 3, with the GCC Runtime Library
  Exception) and, on Windows, MinGW-w64 (its runtime licences are
  permissive).
