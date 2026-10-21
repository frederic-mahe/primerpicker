# Probes of PRIMER.exe, evening of 2026-10-21: the picking screen

Primer/Create -> From Molecule, after its form, opens MODIFY PRIMER
POSITION (so does Modify -> Position): the molecule, five lines of 60
bases, with the primer drawn as a coloured window, its values on one
line, and keys to move it. Recorded with DOSBox-X, headless, `PAUSE` 3,
at most 4 emulators (no fuzzing campaign was running). Recordings and
the exact keys: `picking/` (`README.txt`). Nothing in the README, the
tools or the earlier recordings was changed; new tools:
`tools/ctrlkeys.S` (a DOS TSR: DOSBox-X's `autotype` cannot hold Ctrl,
so F5 to F10 become Ctrl+Home, Ctrl+End, Ctrl+PgUp, Ctrl+PgDn,
Ctrl+Left, Ctrl+Right), `tools/probe_picking.sh` (`probe_keys.sh` with
that TSR, on a drive D: of its own, and the screens' colours),
`tools/screen_colours.py` (each cell's foreground and background, VGA
colour names), `tools/parse_picking.py`, `tools/picking_table.py`
(settled states: `picking/picking.tsv`), `tools/make_width_molecule.py`
(`picking/TWIDTH.GB`).

Predictions: the Find-next positions and verdicts were compared with
primerpicker's Primer Analysis verdict on every position before the
runs (`find_51`, `find_end`, `len_plus`, `comp_400`); the 3'-end run
rule was predicted from Find Pair's ("Find Pair and runs at the 3'
end") and both models were written down before `trun_*` ran.


## Draft for tests/legacy/README.md

## The picking screen (recorded 2026-10-21)

`picking/`: 21 runs of Primer/Create -> From Molecule on
`pairs/TEST600.GB`, `search/TRUN.GB` and `picking/TWIDTH.GB` (written for
these probes, an input), built-in settings; `picking/picking.tsv`: the
71 settled states of the screen (`tools/picking_table.py`); `README.txt`
lists the keys.

**Screen.** "MODIFY PRIMER POSITION", "Primer length: N   Molecule:
NAME", five lines of the molecule (60 bases in blocks of 10, after the
line's first position), "MEETS CRITERIA SELECTED: Yes|No" (with
"*** Complementary Strand ***" on that line for the complement strand,
whose bases are shown complemented), the value line, the keys
"C=Criteria E=Extn(5') I=Inspect M=Mutagenz S=Save+Exit", "Move: (the
arrows) Find next: Ctrl (left), Ctrl (right) Length: + -". Two cells
("--") above the line mark the window's left end. Lines scroll so that
the window stays in view (near the end of TEST600: lines 421 to 541).

**Value line**, for example
`Pos: 51     %GC: 50    Tm°C: 66    Hairpin: OK  Dimers: OK  Runs: !`:
"Pos: " and the position in 7 columns; "%GC:" and "Tm°C:" with the value
right-aligned in 3 columns (`%GC:  5!`, `%GC:100!`, `Tm°C:121!`; a Tm
that is not positive is `Tm°C: ? !`), then `!` when the value is outside
the criteria range or a space, then 3 spaces; "Hairpin: ", "Dimers: "
(3' dimers or adjacent bases) and "Runs: " (runs of any base, or G/C at
the 3' end) followed by `OK` or `!`, padded to 4 columns. **Unlike the
Primer Analysis screen, the run at the 3' end counts** (as in Find Pair,
"Find Pair and runs at the 3' end"): the primer at 450 of TRUN.GB
(GCAGTCAGGCTCATCGCAAA), which meets every Primer Analysis criterion,
shows `Runs: !` and "MEETS CRITERIA SELECTED: No", and Find next skips
it; so do TEST600 at 36 with 21 bases (...CAAA) and AAA (TWIDTH at 41,
3 bases). Predicted before `trun_*` ran (both models written down).

**Colour.** The window's background is green (VGA colour 2) when MEETS
CRITERIA SELECTED is Yes, light grey (7) otherwise, black text in both
(all 71 states); the screen is white on blue; the fields of the value
line keep their colours (black on cyan, the last on light grey)
whatever the verdict. One screen cell (row 12, column 40) keeps a colour
of its own (brown on blue, magenta on green, black on light grey) in
every run; it is not part of the window. With 99 or 100 bases the
window is drawn on 71 bases only (`len_99`).

**Keys** (both strands; the position is the primer's 5' end, as
elsewhere, and grows to the right on both strands):
- Right, Left: position + 1, - 1; Down, Up: + 60, - 60 (one line). A
  move that would take the primer off the molecule is ignored, without
  message (Left at 1, Right at 581 with 20 bases in 600, Up at 51, Down
  at 577; on the complement strand the primer must keep `length` bases
  up to its position).
- +, -: length + 1, - 1 (length 1 is allowed: GC 0, Tm `?`). + that would
  extend past the molecule's end shows "Primer Too Long - Extends
  outside molecule"; - at 1 base shows "Primer must have length!"
  (boxes, "Press any key to continue..."); the primer is kept. The form
  takes 2 digits for the length (99); + then gives 100.
- Ctrl+Right, Ctrl+Left: the nearest position beyond the current one, in
  that direction, where a primer of the same length and strand meets
  every criterion (with the 3'-end run rule above): TEST600 from 51: 71,
  72; back from 72: 71, 36, 35; from 530: 559, then "No Primers found"
  (a box; the primer stays); from 2 backwards: "No Primers found"; with
  22 bases from 36: 71; complement strand from 460(C): 516(C), back from
  516(C): 457(C); TRUN.GB from 440: 449, 471 (450 and 451 end with AAA),
  back: 449; from 40: 42, 46.
- PgUp, PgDn, Home, End, Ctrl+Home, Ctrl+End, Ctrl+PgUp, Ctrl+PgDn: no
  effect.
- S (Save+Exit) **and Esc**: the primer as moved joins the Primer List
  (TEST600 52, `!`; 457(C), `OK`), and the main menu returns.
- Modify -> Position opens the same screen on the active primer
  (`modify_position`: made at 51, saved, moved to 52 and found 71, saved
  again: the list then holds 51 and 71).

Not probed: C, E, I and M on this screen; a primer with a 5' extension;
non-default criteria (the colour follows MEETS CRITERIA SELECTED, which
uses them elsewhere); the "--" mark's glyphs (two cells decoded as
`-`); circular molecules.


## Find Pair's "Press key to interrupt search"

Recordings: `find_pair_interrupt/` (`README.txt`). Three runs typed a
space a few seconds after the last Enter of the form
(`tools/find_pair_limits.sh` with `EXTRA`):
- TEST600, targets 300 and 310, offsets 290, built-in settings: the
  search ended (21 and 34 possible primers, 332 pairs, as primerpicker)
  within a few video frames, before the key; the box then reads
  "<Press key to display possible pairs>", and the space showed the
  list.
- TEST600, targets 150 and 400, offsets 100, `config/PD2_CRITERIA.CFG`
  (the search that reached "Pairs: 1000 possible" and still showed its
  list, "Limits on primers and pairs"): the same, in about 12 frames.
- A random 30,000-base molecule (to slow the search down): keys were
  lost while its map was drawn, twice (not kept).
With the "pairs capped" runs of 2026-10-19 (the program hangs at 1000
pairs and ignores keys), no recorded search could be interrupted:
PRIMER.exe's searches end at once or never. What it shows after an
interrupt is still unknown, and the GUI keeps the behaviour decided on
2026-10-21 ("Search interrupted", no list).


## Mutagenize and the run at the 3' end

`mutagenize_3prime_run/`: Mutagenize on TRUN.GB at 450 (ends with AAA,
meets every Primer Analysis criterion): the ORIGINAL and MUTATED lines
read `Runs: !`, as the picking screen. The six screens of
`mutagenize/mutagenize.tsv` do not tell the two rules apart (no primer
there ends with a run of 3). Predicted from the picking screen.


## Modify -> Position on other primers

Not recorded: a Direct Entry run lost keys ("Not a valid base"). The
executable holds "Directly entered primers cannot be moved" next to
"Mutated primers cannot be moved" (recorded, "Mutated primers"); the GUI
uses these two messages.


## Modify -> Length

Recordings: `modify_length/` (`README.txt`). Modify -> Length opens the
CHANGE PRIMER LENGTH form ("Primer length (incl 5' ext): 20", "Type
over number; Enter To Complete", Esc cancels). With the primer at 51
of TEST600 saved as entry 1 of the Primer List (active):
- 25, Enter: entry 1 becomes `51  25  CTCAAAGGCTCAGACATCAGGGCTC`, in
  place (no new entry); the 5' end stays, the 3' end moves;
- 1: entry 1 becomes `51  1  C`;
- 0: a box "Primer must have length!" (as - on the picking screen);
- 575, 30 bases: a box "Primer Too Long - Extends outside molecule";
  the entry keeps 20 bases;
- 400(C), 25: `400(C)  25  TGTATCCAGGTAATGCAAGGATACC`: on the
  complement strand too, the 5' end (400) stays.
The field takes 2 digits (99 at most; see `picking/len_99`). Not
probed: a primer with a 5' extension (the label says the length
includes it), a Direct Entry or mutated primer.

**Draft for tests/legacy/README.md:** Modify -> Length sets the active
primer's length (1 to 99 bases, the 5' end kept, on both strands), as +
and - of the picking screen do one base at a time, with the same
messages ("Primer must have length!", "Primer Too Long - Extends outside
molecule"); the Primer List entry of the active primer is changed in
place (`modify_length/`).


## Molecule-Ops -> View/Edit

Recordings: `view_edit/` (`README.txt`), on TEST600.GB.
- View/Edit asks for a molecule (Load File), then a menu: Features,
  Sequence, Annotations, Header, Change molecule, ESCape/abandon.
- Every view has a header: `NAME:        TEST600               600
  BPS DNA             UPDATED  10/21/<6`, then `DESCRIPTION:` (empty).
  The year 2026 is written `<6`: the program prints the year minus 1900
  as two digits, the first one past '9' ('0' + 12).
- Sequence: menu bar Edit, Format, Search, Print, Main Menu; lines of 60
  bases in blocks of 10 after the line's first position (right-aligned
  in 8 columns from column 1, then two spaces); Edit opens EDIT SEQUENCE
  (OVERWRITE; INS toggles insert/overwrite; Esc exits/saves). Format
  lists Single Strand, Double Strand, Translate, Block Size Toggle, AA
  code toggle; no choice could be made in the recordings (arrows and
  initial letters were ignored).
- Search -> FIND SEQUENCE: "Search for:" and the IUPAC codes (N any; M
  A/C, R A/G, W A/T, S C/G, Y C/T, K G/T; B not A, D not C, H not G, V
  not T); the result: `Number of Sites =  N` and the sites on **both
  strands**, by position (the leftmost base, `C` after a complement
  site: TTAAAT gives `8` and `589C`, its reverse complement ATTTAA being
  at 589), each in a field of 10 columns, six a line, comma-separated
  (AAAT: 10 sites on two lines); IUPAC letters match (CATNTA: 5 sites).
  These are the sites of the restriction-site search (true sites),
  without mismatches. Predicted for TTAAAT and CATNTA from that search.
- Features: `-- no features --` (TEST600 has none); keys: Change, INS =
  add feature, DEL = remove.
- Annotations: the lines of the file before the sequence, ORIGIN
  included (LOCUS, DEFINITION, FEATURES, ORIGIN).
- Header: EDIT MOLECULE HEADER: `Name: TEST600`, `Size: 600     sequence
  present`, `Description:` (empty: the DEFINITION line is not used).
- Save/Remove with no molecule loaded: "Molecule List is Empty"; with
  TEST600 loaded: SAVE/REMOVE MOLECULE, the list of loaded molecules,
  "[ S = Save ]    [ R = Remove ]", and "Removing a molecule from the
  list unloads the molecule from program memory (it does not delete a
  disk file).  If the molecule has just been created, save to disk
  before removing the file." What S writes was not recorded.

**Draft for tests/legacy/README.md:** the paragraph above, as a section
"Molecule-Ops -> View/Edit (recorded 2026-10-21)".
