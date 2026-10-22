# Probes of PRIMER.exe, 2026-10-22: saving and creating molecules

Molecule-Ops -> Save/Remove -> S, and Molecule-Ops -> View/Edit ->
Create New. Recorded with DOSBox-X, headless, `PAUSE` 3 to 3.5, at most
4 emulators (two runs froze early with 4; later batches used 2 or 3).
Recordings and exact keys: `molecule_save/` (`README.txt`). Nothing in
the README, the tools or the earlier recordings was changed. The
format was read from the files PRIMER.exe wrote, then checked by
loading hand-written files (`molecule_save/inputs/`), designed so that
competing readings disagree.


## Draft for tests/legacy/README.md

## Saving and creating molecules (recorded 2026-10-22)

`molecule_save/` (`README.txt`): 27 runs; `written/` holds every file
PRIMER.exe wrote, `inputs/` the files written for the probes.

**Save.** Save/Remove -> S on a loaded molecule opens "Save as
filename:" with `C:\NAME.EXT` of the file it was loaded from (the
field shows `?:\` for a moment), "(current directory: C:\)". Enter on
an existing file asks "File already exists!  Overwrite (Y/N)?"; Y
overwrites it **in PRIMER.exe's own format**, whatever the file's
format was (`save_fa`: TFASTA.FA becomes a native file); typing a name
(`ZZ.GB`) writes a new file, in the same format (the extension does
not choose the format). A molecule made with Create New is offered
"Molecule Created. Save a Copy to Disk (Y/N)?", then "Save as
filename:" (empty); the name is written in upper case, without an
extension (the executable holds "Filename must not include an
extension"; `zn` gives `ZN`); Enter on the empty field saves nothing.

**The native format** (DOS line ends, one value a line):
1. the name (from LOCUS or the FASTA header: TEST600, TCIRC; as typed
   for Create New: `zn`);
2. the size (the number of bases, or the Size field of Create New);
3. `0` in every file written (a hand-written `1` gives "Illegal File
   Format", with two record layouts tried: `hand1`, `hand9`, `handA`);
4. the number of GENE features, then 5 lines per gene: `NAME` NUL
   `DESCRIPTION` on one line, start, end, `cw` or `ccw`, `Y` or `N`;
   `ccw` is written when the features table shows `(C)` for the gene
   (below), `cw` otherwise, whatever the file read had (`round_zh6`:
   a gene 2 to 6, `ccw`, `N` is saved `cw`; `round_zh7`: 2 to 6,
   `ccw`, `Y` stays `ccw`);
5. the number of MARKER features, then 3 lines per marker: `NAME` NUL
   `DESCRIPTION`, position, `R` or `L` (`L` when Complement is Y);
6. linear molecules only: `linear`, the left end, the right end, then
   two lines that are `^` and `^` when both ends are empty and empty
   lines when they were typed (`e`, `b`); circular molecules have none
   of these lines (`save_circ`);
7. when there is a sequence: `sequence` (DNA) or `sequenceAA`
   (protein), then the whole sequence on one line, in upper case;
8. when the description is not empty: `description`, then the text.
   PRIMER.exe reads the description before or after the sequence
   (`hand3`, `hand7`), and writes it after (`round_zh7`: ZH7, with the
   description before the sequence, saved over itself, has it after).
A protein made with Create New (P) is written as DNA would be
(`create_prot`: `0`, `linear`, no sequence): only `sequenceAA` marks a
protein (`hand4`: "4 AA PROTEIN"). The description of a GenBank file
(DEFINITION) or of a FASTA header is not kept (empty description).

**Features.** REGION features made with Create New are listed, but
not saved (`feat_region2`: 0 genes, 0 markers). FEATURES lines of a
GenBank file are not read (`tfeat_view`: "-- no features --"). In the
features table, MARKER rows show `»` (R) or `«` (L) in the End column,
and `(C)` for L. A gene's `(C)` follows neither `cw`/`ccw` nor the
flag alone: it is shown when the flag is `Y` and start < end, or `N`
and start > end (`hand6`, `hand7`, `hand8`: 6 genes); a gene typed
from 9 to 5 with Create New is written `cw`, `Y` and shown without
`(C)` (`feat_gback`). How a gene gets Complement in the Create New
form was not found (Down after End base did not reach it).

**Reading.** PRIMER.exe reads the files it writes (`reload_header`,
`reload_primer`: name TEST600, 600 bases, a primer can be made from
it; `reload_zn`: the gene). The header view of a native file shows
`Name: TEST600      :` and, for circular files, `CIRCULAR` after
`BPS DNA`; a description is shown (`hand3`, `hand7`).

**Create New** (`create_open`): CREATE NEW MOLECULE, "Name:", "Size:",
"DNA/Protein (D/P): D", "Circular (Y/N)  : N", "Left End:", "Right
End:" (both hidden for proteins), "Description:"; then "Enter Features
(Y/N)?" (the features table: INS adds a feature, form INSERT/EDIT
FEATURE: name, description, type G/M/R, start, end — empty for
markers —, complement), "Enter Sequence (Y/N)?", a summary ("Number of
Features Entered = N", "Sequence Entered:  No") and the save prompt.
Answering Y to "Enter Sequence" stopped the program in every run
(`create_seq4`: the screen stays on "Sequence Entered:  No" and keys
are ignored): a sequence could not be typed.

Not probed: Remove (R) after a save; a description and a sequence
written by Create New (the editor hangs); a molecule with a sequence typed
in (the editor hangs); left and right ends named after enzymes; the
meaning of the two `^` lines; how a created gene gets Complement;
whether `cw`/`ccw` matters on circular molecules.
