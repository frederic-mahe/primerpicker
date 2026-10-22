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
