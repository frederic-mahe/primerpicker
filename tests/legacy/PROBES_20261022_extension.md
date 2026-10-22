# Probes of PRIMER.exe, 2026-10-22: primers with a 5' extension

Open questions of the proposal (§6.3f, §6.4): Modify -> Length and the
picking screen on a primer with a 5' extension; uniqueness, Primer Pair
and Multi Primer with extensions. Recorded with DOSBox-X, headless,
`PAUSE` 3.5, 1 or 2 emulators. Recordings and keys: `extension/`
(`README.txt`). Nothing in the README, the tools or the earlier
recordings was changed; new tool: `tools/picking_table_extension.py`
(`picking_table.py` for a window shorter than the primer).

Predictions: the value lines were compared with `primerpicker analyze
--extension GAATTC CTCAAAGGCTCAGA` (GC 45, Tm 47, runs 3: `!`) before
the comparison.


## Draft for tests/legacy/README.md

## Primers with a 5' extension (recorded 2026-10-22)

`extension/` (`README.txt`): a primer made with From Molecule at 51 of
TEST600, length 20 ("Primer length (incl 5' ext)"), extension GAATTC:
`GAATTC-CTCAAAGGCTCAGA`, 14 bases of the molecule.

**Picking screen** (`ext_picking3`, `picking_extension.tsv`, 5 settled
states). The header reads `Primer length: 20   Molecule: TEST600
5'extension: GAATTC`: the length counts the extension. The value line
gives the values of the whole primer, the Tm without the extension, as
the Primer Analysis screen (`Pos: 51     %GC: 45!   Tm°C: 47!
Hairpin: OK  Dimers: OK  Runs: !`). The window covers the molecule's
bases only (14 bases from the position). + gives `Primer length: 21`
(15 bases of the molecule: 52, `%GC: 42!   Tm°C: 52!`); Ctrl+Right
from 52 finds 75 (`%GC: 50    Tm°C: 56`, every criterion met, the
window green): Find next uses the values of the whole primer.

**Modify -> Length** (`ext_length2`): the form shows `Primer length
(incl 5' ext):  20`; 25 changes entry 1 of the Primer List to
`51  25  GAATTC-CTCAAAGGCTCAGACATCA` (the extension and 19 bases of the
molecule).

**Primer Pair** (`ext_pair`: A as above, B at 400(C), 20 bases): the
PRIMER PAIR ANALYSIS screen gives primer A's values as the Primer
Analysis screen does (length 20, `%GC 45!`, `Tm 47!`: the whole primer,
the Tm without the extension; `19°C difference`), and the A:B dimers
of the whole primers (3' end 3!, any 3). The amplified product is that
of the primers without extension (51 to 400: 350 bases, Tm 77, as
`pairs/pairs.tsv` for 51 and 400C), and the annealing temperature uses
primer A's Tm without extension: floor(0.3 x (47 - 21.6) + 0.7 x 77 -
14.9) = 46.

**Multi Primer** (`ext_multi`: M13EXT.PD2, GAATTC +
GTAAAACGACGGCCAGT, and AACAGCTATGACCATG): the same rule (A: length 23,
`%GC 47!`, `Tm 64`, runs 4!, 3' dimers A:A 2, any 6); `primerpicker
multi` on the whole sequence gives the same values except the Tm (77).

**Uniqueness** (`ext_unique`): "Primer is unique" (as `primerpicker
homology` for 51 without extension): this run does not tell whether
the extension is searched.

**Reports** (`pair_ext_report.PRN`, `multi_ext_report.PRN`: Print to
Disk of the Primer Pair and Multi Primer screens of M13EXT.PD2 and
AACAGCTATGACCATG): the sequence is written with a dash after the
extension (`Sequence:     5'-GAATTC-GTAAAACGACGGCCAGT-3'`), the length
and composition are those of the whole primer (`Length:  23;`,
`A:8  C:5  G:6  T:4`), the Tm is without the extension (64), and the
duplexes of page 2 are drawn on the whole primers, without dash. No
other line mentions the extension.
