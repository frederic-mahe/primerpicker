> **Note:** primerpicker is developed with the help of Anthropic's
> Claude (Opus 5.5).

# primerpicker

graphical picker for DNA primers

primerpicker designs and checks PCR and sequencing primers. It is
inspired by *Primer Designer* v2.0 (Scientific & Educational Software,
1990-1991), a graphical tool for DOS that some laboratories still run in
an emulator. Version 1 reproduces its features, default values and
results; later versions will add modern methods.

- implemented in modern C++ (C++20), with a core library that has no
  user interface and is unit-tested;
- aiming for high-quality programming standards: test-driven
  development, warning-free builds with GCC and clang, sanitizers,
  static analysis (clang-tidy, cppcheck), fuzzing of every reader of
  untrusted input;
- a graphical user interface *and* a command-line interface (the
  graphical interface, in Qt 6, is under development; the
  command-line interface is usable now);
- runs on Windows, macOS and GNU/Linux.

Every result is checked against reference values: screens recorded
from Primer Designer itself (`tests/legacy/`, with the tools that
recorded them), published tables, and reference implementations
(Primer3, REBASE).


## What it can do

Command-line interface (`primerpicker --help`):

- **`analyze`**: one primer, with an optional 5' extension: length,
  base composition, GC %, melting temperature, hairpins, self-dimers
  (matches at the 3' end, adjacent homologous bases), runs of bases
  and of G/C at the 3' end, each checked against the criteria;
- **`pair`**: two primers: balance (length, GC % and Tm differences),
  cautions, and the A:A, B:B and A:B dimers;
- **`multi`**: two to four primers: balance, cautions, and the dimer
  matrices of the Multi Primer Analysis;
- **`find-pair`**: pairs of primers around two targets of a molecule,
  ranked with five weighting factors (GC and Tm mismatch, 3' dimers,
  any dimers, false priming), with chosen search windows;
- **`find-enzymes`**: restriction sites on both strands of a linear or
  circular molecule, exact or with one mismatch, for the commercially
  available enzymes of REBASE (old enzyme names are accepted);
- **`--config`**: the criteria, tolerances, search windows and weights
  of a Primer Designer settings file (`pd2.cfg`).

Also in the core library, for the graphical interface: Match Active
Primer (a partner for a given primer), uniqueness and false priming
on the template, PCR product Tm and annealing temperature, and the
molecule (GenBank, FASTA, plain sequence) and primer (`.PD2`) files
of Primer Designer.


## Molecule size

primerpicker has no fixed limit on the length of a molecule (Primer
Designer refused a sequence text of 64 KB or more: about 51,000 bases
in a GenBank file, `tests/legacy/README.md`). Positions are
`std::size_t` throughout (64-bit on the supported platforms); the
largest molecule tested has 4.3 billion bases,
more than 2^32 (`tests/large/`, an opt-in test). In practice the limit
is memory and time. Measured on an Intel Core Ultra 9 285 (release
build, one thread):

| Molecule | Reading + `find-enzymes` | Memory | `find-pair` | `find-pair`, 99 + 99 primers |
|---|---|---|---|---|
| 1 million bases | 0.01 s | 8 MB | 3 s | 28 s |
| 10 million bases | 0.14 s | 40 MB | 29 s | 279 s |
| 100 million bases | 1.3 s | 330 MB | 217 s | (not measured) |
| 4.3 billion bases | 54 s | 12.7 GB | (not measured) | (not measured) |

Memory is about 3 bytes per base while a molecule is read. `find-pair`
searches the whole molecule for false-priming sites of each possible
primer: its time grows with the molecule and with the number of
possible primers (random molecules, targets 300 and 900, built-in
criteria; last column: relaxed criteria, so that 99 primers qualify
on each side, the most a search accepts).


## Example

```console
$ primerpicker analyze GTAAAACGACGGCCAGT
Sequence:      5'-GTAAAACGACGGCCAGT-3'
Length:        17
5' extension:  --none--
Composition:   A:6  C:4  G:5  T:2

Criteria         Setting                                  Actual  Meets
% GC             min 50, max 60                               52  YES
Tm °C            min 55, max 80                               64  YES
No hairpins      energy cutoff 0.0 kcal                       --  YES
No 3' dimers     reject >= 3 matches at 3' end                 2  YES
No dimers        reject >= 7 adjacent homologous bases         4  YES
No runs          reject >= 3 base runs                         4  NO
No 3' GC runs    reject >= 3 G or C at 3' end                  0  YES

Hairpin:       no hairpins found
Meets criteria: NO
```


## Building

Requirements: CMake 3.28 or later, Ninja, and a C++20 compiler (GCC 13
or later, clang 18 or later, or Apple clang). Catch2 3 is used when it
is installed, and downloaded otherwise.

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

The `release` preset builds with GCC; the program is
`build/release/src/cli/primerpicker`. Other presets
(`CMakePresets.json`) build with sanitizers (`debug`), with clang and
libc++ (`clang-debug`), for Windows with MinGW-w64 (`mingw-release`),
and the fuzzing harnesses (`fuzz`).

Pre-compiled executables will be attached to the releases.


## Citation

See `CITATION.cff` (GitHub shows it as "Cite this repository").


## License

primerpicker is free software, released under the GNU General Public
License, version 3 or later (`LICENSE`).

primerpicker is not affiliated with Scientific & Educational Software;
Primer Designer is their product, and is not part of this repository.

The restriction enzyme data come from REBASE (version 610, files in
`data/rebase/`): Roberts RJ, Vincze T, Posfai J, Macelis D (2010)
REBASE - a database for DNA restriction and modification: enzymes,
genes and genomes. *Nucleic Acids Research* 38: D234-D236.
[doi:10.1093/nar/gkp874](https://doi.org/10.1093/nar/gkp874)
