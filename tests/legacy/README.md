# Recorded behaviour of the legacy program

Reference observations of *Primer Designer v2.0* (Scientific &
Educational Software, 1990-1991), `PRIMER.exe`. The executable itself
is proprietary and is not part of this repository.

Recordings are made with DOSBox-X 2026.01.02, running headless:
keystrokes are injected with the built-in `AUTOTYPE` command, and
the screen is recorded with `dx-capture /V` (ZMBV video). The program
runs in 80x25 text mode (720x400 pixels).

- `tools/zmbv.py`: decodes a ZMBV AVI file into PNG images of the
  distinct frames (needs Python 3 and Pillow).
- `tools/system_defaults.conf`: DOSBox-X configuration that records
  the System Defaults screens (usage in the file header).
- `tools/screen_text.py`: decodes the PNG frames into text, exactly
  (each 9x16 cell is matched against the VGA font in
  `tools/vga_glyphs.png`; `tools/make_vga_glyphs.sh` remakes it).
- `tools/tm_direct_entry.sh`: types one primer with Primer/Create ->
  Direct Entry and prints the result line (usage in the file header).
- `tools/fit_tm.py`: checks the primer Tm model below against
  `tm/direct_entry_tm.tsv`.
- `tools/make_molecule.py`: writes a random test molecule as a
  GenBank file.
- `tools/pair_product.sh`: makes two primers from a molecule, and
  records the Primer Pair and amplified product screens.
- `tools/fit_pairs.py`: checks the product Tm and annealing
  temperature models below against `pairs/pairs.tsv`.
- `tools/homology_screen.sh`: makes a primer from a molecule and
  records its Uniqueness or False Priming screen;
  `tools/parse_homology.py` turns the screen into table rows.
- `tools/find_pair.sh`: records the Find Primer Pair list for two
  targets on a molecule, with optional weights.
- `tools/multi_primer.sh`: loads up to four primers from `.PD2`
  files (written by `tools/make_pd2.py`) and records the Multi Primer
  Analysis screen; `tools/parse_multi.py` turns it into tables.
- `tools/find_pair_limits.sh`: runs Find Primer Pair with chosen
  search windows and settings file, and prints the counts and
  messages (limits on primers and pairs).
- `tools/pair_direct.sh`: types two primers with Direct Entry and
  records their Primer Pair Analysis screen; `tools/parse_pair.py`
  turns pair screens (from it or from `pair_product.sh`) into a table.
- `tools/inspect_primer.sh`: types one primer, and records the Primer
  Analysis (Inspect), Hairpins and Dimers screens;
  `tools/parse_inspect.py` turns the recordings into a table.


## Built-in defaults (recorded 2026-10-17)

Fresh start, no `pd2.cfg` file present. Screens: Modify -> System
Defaults.

| Screen | Setting | Default |
|---|---|---|
| `screens/system_defaults.png` | Primer length (incl. 5' extension) | 20 bases |
| | False priming and uniqueness search: exclude Tm <= 0 °C | Y |
| | Sort homologies by Tm | Y |
| | Match primer GC% to within ± | 5 % |
| | Match primer Tm to within ± | 10 °C |
| | Match Active Primer search range | start -50, end +50 bases from target |
| | Find Pair search range | normal start -100, complement end +100 bases from target |
| `screens/criteria_defaults.png` | GC% range | 50-60 % |
| | Tm range | 55-80 °C |
| | Hairpins: energy cutoff | 0.0 kcal |
| | Dimers: reject if >= | 3 matches at 3' end |
| | Dimers: reject if >= | 7 adjacent homologous bases |
| | Runs: reject runs of >= | 3 bases (any) |
| | Runs: reject runs of >= | 3 G or C at 3' end |
| `screens/weights_defaults.png` | Find-Pair weighting factors (0 = ignore, 1-5): GC mismatch, Tm mismatch, 3' dimers, any dimers, false priming | 1 each |


## Primer Tm (recorded 2026-10-17)

`tm/direct_entry_tm.tsv`: 186 primers typed with Primer/Create ->
Direct Entry (built-in defaults), and the values shown on the result
line (`Length: 17  %GC: 52    Tm°C: 64 ...`). A `!` flag means that
the value is outside the criteria range. The primers are homopolymers
and dinucleotide repeats (12 to 30 bases), random primers (8 to 64
typed bases, 20 to 80 % GC), primers chosen to narrow the fit below,
primers whose Tm is close to 0 °C, five primers with a 5' extension,
and two with non-ACGT letters.

All 186 values follow this model (`tools/fit_tm.py`):

    Tm = floor(1000 * dH / (dS + b) - 273)    displayed as "?" when Tm <= 0

- dH (kcal/mol) and dS (cal/K/mol) are sums over the nearest-neighbour
  pairs of the primer, with the parameters of Breslauer et al. 1986
  (PNAS 83:3746, Table 2). No initiation term, no terminal or
  symmetry correction (self-complementary repeats such as `ATAT...`
  and `CGCG...` follow the same formula).
- b = 60.698 to 60.700 cal/K/mol, constant (no length, salt or
  concentration term was needed). Its origin is not known. One
  combination that falls in the interval: Breslauer's initiation free
  energy (5.0 kcal/mol) divided by 298.2 K, plus R ln(1/C) with
  R = 1.987 cal/K/mol and C = 250 pM: 16.767 + 43.932 = 60.699. This
  is a guess, not a finding.
- The result is truncated (floor), and the offset is 273, not 273.15:
  rounding, or 273.15, contradict several records.
- %GC is truncated as well (9 G+C out of 17: 52.9 % shown as 52).
- The 5' extension does not change the Tm, but it is included in the
  displayed length and %GC.
- The primer field holds 50 bases; further keys are ignored.
- Letters other than A, C, G and T (N, R, U) are not inserted, and
  the next key is lost as well; `U` garbles the entry field. Not
  investigated further.

Nearest-neighbour parameters were not found in the executable (as
float, double or 16-bit integers). The executable does contain these
doubles, next to each other at file offset 0x429d3: 21.6, 0.3, 41,
675, 59.9, 0.7, 14.9. They look like a product Tm of
59.9 + 41 x (GC fraction) - 675 / length (that is, 81.5 + 16.6
log10(0.05) + 0.41 %GC - 675 / length) and the annealing temperature
of Rychlik et al. 1990 (0.3 x Tm primer + 0.7 x Tm product - 14.9).
Confirmed with Primer Pair analysis: see "Product Tm and annealing
temperature" below.

Recording notes:
- Keys are typed every 0.3 s. On a busy machine (about 20 runs in
  parallel), the program misses keys (probably the 15-key BIOS
  buffer); slower typing (`PAUSE=1`) fixes it. Each record was
  checked: the displayed length equals the typed length (except
  where the 50-base limit or non-ACGT letters apply).
- DOSBox-X resolves a relative `mount c c` against the directory of
  the configuration file, not the current directory.


### Primer3 equivalent (2026-10-18)

Primer3's `oligotm` (2.6.1, Breslauer mode: `-tp 0 -sc 0 -dv 0 -n 0`)
gives the same Tm, to its 6 printed decimals, for all 159 recorded
primers it accepts (A, C, G and T only, up to 36 bases), with:

- DNA 0.0496293659554817 nM: Primer3 computes dH / (dS - 10.8 +
  R ln(C/4)) with R = 1.987, so b = 10.8 + R ln(4/C) = 60.699 for
  C = 49.6 pM. A primer concentration of about 50 pM is therefore one
  possible origin of b;
- monovalent cations 1021.0244559195794 mM: this has no physical
  meaning, it makes Primer3's salt term 16.6 log10([K+]) equal to
  +0.15 so that its -273.15 becomes Primer Designer's -273.

The values are in `../data/primer3_legacy_tm.tsv` (made by
`../data/primer3_legacy_tm.sh`) and tested in
`tests/unit/test_legacy_tm_primer3.cpp`.

## Product Tm and annealing temperature (recorded 2026-10-17)

`pairs/pairs.tsv`: 28 primer pairs on two generated molecules, made
with Primer/Create -> from Molecule (primer A at a position, 20 bases;
primer B on the complement strand), then Analyze -> Primer Pair ->
A=AnalProd. Molecules (GenBank files, which the program loads):
- `pairs/TEST600.GB`: 600 random bases (Python `random.Random(600)`,
  `choice('acgt')`).
- `pairs/TGRAD.GB` (LOCUS name GRAD): 1,200 bases in eight segments
  of 150 with GC fractions 0.25, 0.4, 0.55, 0.7, 0.55, 0.4, 0.25,
  0.8 (`tools/make_molecule.py GRAD 1201
  0.25,0.4,0.55,0.7,0.55,0.4,0.25,0.8 150`).

The product spans primer A's 5' end to primer B's 5' end (positions
51 and 400C give 350 bases). All 28 records follow
(`tools/fit_pairs.py`):

    product Tm = floor(59.9 + 41 x GC fraction - 675 / length)
    Ta = floor(0.3 x (lower primer Tm - 21.6) + 0.7 x product Tm - 14.9)

- The constants are the doubles found in the executable (see "Primer
  Tm" above): 59.9 = 81.5 + 16.6 log10(0.05), and 21.6 =
  -16.6 log10(0.05), so both formulas assume 50 mM monovalent salt.
- Ta uses the *displayed* (truncated) Tm values, of the primer with
  the lower Tm (A or B, both occur in the data), and of the product.
  Using unrounded values, or another primer choice, contradicts
  several records.
- The GC fraction is not truncated to a whole percentage (that
  contradicts 7 records).
- The primer Tm values on the pair screen are those of Direct Entry.

Not yet probed: products with a 5' extension, Tm values at or below
0, and primer lengths other than 20.


## Primer criteria (recorded 2026-10-17)

`inspect/inspect.tsv`: 111 primers typed with Direct Entry, then
Inspect (Primer Analysis screen, "Actual" column), H (Hairpins screen)
and D (Dimers screen), with built-in defaults. The primers are those
of `tm/direct_entry_tm.tsv` (first batch) whose recording completed.

Counts, exact for every record (`tests/unit/test_legacy_criteria.cpp`
checks them all):
- **Runs**: the longest run of one base, not counting the run at the
  3' end (a homopolymer gives 0; `CAACTGCTATTT` gives 2).
- **3' GC runs**: the number of G or C at the 3' end.
- **Dimers, adjacent homologous bases**: the longest run of
  consecutive complementary pairs between the primer and a copy of
  itself, in any antiparallel alignment (the fully overlapping one
  included).
- **3' dimers, matches at 3' end**: in antiparallel alignments of the
  primer with a copy of itself where the copy extends at least one
  base beyond the primer's 3' end (a template for extension), the
  largest number of consecutive complementary pairs counted from the
  primer's 3' end. Without that overhang condition, 19 records
  disagree (for instance, self-complementary 30-mers give 28, not 30).

**Verdict** ("Meets Criteria" column, checked on 184 recorded Primer
Analysis screens): GC% and Tm pass inside their ranges, limits
included (YES is shown at 50 and 60 %, and at 55 and 80 °C); the
four counts fail when they reach the "Reject >=" value; the hairpin
criterion fails when a hairpin is found. "MEETS CRITERIA SELECTED" is
YES exactly when all seven criteria are met. Where Direct Entry shows
the Tm as "?", the Primer Analysis screen shows 0.

**Hairpins** (Hairpins screen): an energy in kcal with one decimal
(`-.4` for -0.4), and a drawing of one hairpin: the 5' part, pairing
marks (`:`), and the 3' part written 3' to 5', right-aligned, with
`)` for the loop; or "No Hairpins Found". When the 3' part is long,
its left end is cut off the screen; for an odd loop, the middle loop
base is drawn on both lines.

`inspect/hairpin_design.tsv` adds hairpins designed to vary one thing
at a time: stem length, loop length (3 to 12) and composition, single
nearest-neighbour stems, and flanking bases.

Hairpin model (every recorded hairpin, energy and "No Hairpins
Found", is reproduced by `tests/unit/test_legacy_hairpin.cpp`):
- energy = loop penalty - sum of stacking stabilities along the stem,
  in kcal, shown with one decimal;
- stacking stabilities: Freier et al. 1986 (PNAS 83:9373), RNA
  nearest neighbours at 37 °C, read for DNA (T for U): AA 0.9, AT
  0.9, TA 1.1, CT 1.7, CA 1.8, AC 2.1, GA 2.3, CG 2.0, GG 2.9, GC 3.4
  (and their equivalents on the other strand). Not Breslauer 1986,
  which the primer Tm uses;
- loop penalty: Freier 1986 hairpin loops, by length: 3: 7.4, 4: 5.9,
  5: 4.4, 6: 4.3, 7: 4.1, 8: 4.1, 9: 4.2, 10: 4.3, 12: 4.9, 14: 5.6,
  16: 6.1, 18: 6.7, 20: 7.1, 30: 8.9. Lengths in between are
  interpolated linearly (11: 4.6, 13: 5.25, 23: 7.64), and lengths
  beyond 30 extrapolated with the slope from 20 to 30 (0.18 per base;
  recorded up to 34: a loop of 35 next to the strongest 5-pair stem
  gives no hairpin). The loop sequence and bases outside the hairpin
  do not matter;
- the stem may contain G.T pairs (as RNA G.U), but never at either
  end: a stem starts and ends with Watson-Crick pairs. Steps with a
  G.T pair have their own values; the recordings determine only the
  sum of the two steps around each G.T pair, so the ten step values
  in `core/thermo/freier_1986.hpp` are one solution of that system;
- stems are zipped inwards as long as the next inner pair is a
  Watson-Crick pair and the loop keeps at least 3 bases; only these
  maximal stems are scored (`ATAT...` repeats get a loop of 4, not
  6);
- the most stable hairpin is shown; "No Hairpins Found" when none has
  a negative energy;
- the energy is shown truncated toward zero to tenths: a loop of 23
  gives -2.16, shown -2.1 (all values are exact in hundredths of
  kcal, so the core computes in integer hundredths);
- stems do not extend through a tandem T.G/G.T step (no recorded
  stem contains one).

`inspect/inspect_2.tsv` adds the primers of the first batch whose
recording had failed, loop lengths 13 to 40 next to the same stem,
hairpins with one or two G.T pairs inside 6- or 7-pair stems, and
G.T pairs offered at the outer end of a stem or next to the loop
(they are not used).


## Molecule file formats (recorded 2026-10-18)

Files in `molecules/` hold the sequence of `pairs/TEST600.GB` in other
formats (DOS line ends). Each was loaded with Primer/Create -> from
Molecule -> Load File, then a primer was made at position 51 (20
bases), which shows the first sequence lines and the primer values.

| File | Format | Result |
|---|---|---|
| `pairs/TEST600.GB` | GenBank (`LOCUS`, `ORIGIN`, `//`) | loads; name from `LOCUS` (TEST600) |
| `molecules/TFASTA.FA` | FASTA, uppercase, header `>TEST600 random test molecule` | loads; name TEST600 |
| `molecules/TFA2.FA` | FASTA, lowercase, header `>VERYLONGMOLECULENAME_2026 with a description` | loads; name VERYLONGMOLECULENAME_2026 (first word, not truncated) |
| `molecules/TFA3.FA` | FASTA, uppercase, header `>name with spaces` | loads; name `name` |
| `molecules/TPLAIN.SEQ` | plain sequence, no header | loads; name TPLAIN.SEQ (the file name) |
| `molecules/TEST600.EMB`, `TMIN.EMB`, `TNONUM.EMB` | EMBL (`ID`, `SQ`, `//`), with or without other lines or position numbers | "Illegal File Format" |

Every file that loads gives the same sequence (first line
`GGTCATTTTA AATGGACAAC ...`) and the same primer at position 51
(%GC 50, Tm 66). The executable contains the strings `ID   `,
`FT   ` and `SQ   `: EMBL may be written (Molecule-Ops -> Save) or
read in a form not found yet.


## Molecule size (recorded 2026-10-19)

`molecules/molecule_size.tsv`: molecule files of increasing size
loaded with Primer/Create -> from Molecule -> Load File
(`tools/load_molecule.sh`). The files are not kept: each row gives the
command that writes it again (`tools/make_molecule.py`, one segment of
GC 0.5, seed 7; `tools/make_sized_fasta.py`; or, for two FASTA rows,
the `make_molecule.py` sequence written as FASTA, 60 bases per line,
DOS line ends), its size, and the size of its sequence text (the
lines after `ORIGIN`, or after the FASTA header line, line ends,
spaces and numbers included).

The limit is on the sequence text, not on the number of bases or on
the file size:
- 65,536 bytes of sequence text or more: "Sequence too large" (and
  "Press any key to continue..."), the molecule is not loaded;
- just below, the load may fail for lack of memory: 65,534 and 65,535
  bytes give "Insufficient Memory Error" (DOSBox-X default memory
  settings), 65,450 bytes and less load;
- so the largest molecule depends on the format: about 51,000 bases
  in a GenBank file (about 1.28 bytes of text per base: numbers and
  spaces), about 63,000 bases in a FASTA file with DOS line ends.
primerpicker does not reproduce this limit (decided 2026-10-19:
documented, not reproduced; its own limits are in the top-level
README, "Molecule size").


## Settings file `pd2.cfg` (recorded 2026-10-18)

The program writes `PD2.CFG` (55 bytes) when it quits after System
Defaults were changed ("System Defaults Changed But Not Saved! SAVE
CURRENT DEFAULTS (Y/N) ?", answered Y).
- `config/PD2_LENGTH25.CFG`: built-in defaults, primer length set to
  25.
- `config/PD2_CRITERIA.CFG`: primer length 20; criteria changed to GC
  45-65 %, Tm 50-75 °C, hairpin cutoff -1.5 kcal, 3' dimers 4,
  adjacent bases 8, runs 4 (the 3' GC run stayed at 3: the key was
  lost); weights changed, but keys were lost there too.
- `config/PD2_WEIGHTS12345.CFG` and `config/PD2_WEIGHTS35142.CFG`
  (recorded 2026-10-19 with `tools/save_weights.sh`): built-in
  defaults, with the weights of the form (GC mismatch, Tm mismatch,
  3' dimers, any dimers, false priming) set to 1 2 3 4 5 and to 3 5 1
  4 2. Bytes 42-46 hold 1 2 5 3 4 and 3 5 2 1 4: the false-priming
  weight is stored third.
- `config/PD2_RUNS42.CFG` (recorded 2026-10-19 with
  `tools/save_run_thresholds.sh`): built-in defaults, with "Reject
  runs of >= 4 bases (any)" and "Reject runs of >= 2 G or C at 3'end".
  Offset 38 holds 4 and offset 40 holds 2: the order below is
  confirmed. In the Criteria form the dimer fields take two digits
  (Enter moves on) and the run fields one digit (they move on by
  themselves; the last one completes the form).

Layout (little-endian; int16 unless noted), confirmed by the two files:

| Offset | Value | Default |
|---|---|---|
| 0 | unknown (4 bytes, zero in both files) | 0 |
| 4 | primer length (incl. 5' extension) | 20 |
| 6 | exclude Tm <= 0 °C (1 = Y) | 1 |
| 8 | sort homologies by Tm (1 = Y) | 1 |
| 10 | match GC to within ± (%) | 5 |
| 12 | match Tm to within ± (°C) | 10 |
| 14, 16 | Match Active Primer range: start -, end + (bases) | 50, 50 |
| 18, 20 | Find Pair range: normal start -, complement end + | 100, 100 |
| 22, 24 | GC min, max (%) | 50, 60 |
| 26, 28 | Tm min, max (°C) | 55, 80 |
| 30 | hairpin energy cutoff (float32, kcal) | 0.0 |
| 34 | reject 3' dimers if >= | 3 |
| 36 | reject adjacent homologous bases if >= | 7 |
| 38 | reject runs if >= | 3 |
| 40 | reject 3' G/C runs if >= | 3 |
| 42-46 | five weights, one byte each: GC mismatch, Tm mismatch, false priming, 3' dimers, any dimers (recorded 2026-10-19, below) | 1 |
| 47-54 | zero in both files | 0 |

Which of offsets 38 and 40 holds which run threshold was first
inferred from the screen order (only offset 38 changed);
`config/PD2_RUNS42.CFG` confirms it (recorded 2026-10-19). Navigation notes: in the
settings forms the first digit typed clears a field; Enter on the
last field completes the form, and the menu then reopens the same
form by itself (Esc cancels that second pass).


## Primer pair analysis (recorded 2026-10-18)

`pairs/pair_screens.tsv`: the Primer Pair Analysis screen (balance,
cautions, dimer matrix) of 36 pairs: those of `pairs/pairs.tsv`, and
pairs typed with Direct Entry, where B's 3' end is complementary to
A's 3' end over 0 to 6 bases. `pairs/pair_screens_2.tsv` adds 6
pairs, four of them designed to have exactly 7 adjacent homologous
bases (A:B or A:A).

All 36 follow these rules (`tests/unit/test_legacy_pair.cpp`):
- self values (A:A, B:B) and cautions are the single-primer counts
  (see "Primer criteria");
- A:B adjacent homologous bases: the longest run of complementary
  pairs between A and B in any antiparallel alignment;
- A:B 3' end matches: the larger of A's 3'-end count on B and B's
  3'-end count on A (each with the partner extending beyond the 3'
  end, as for self-dimers);
- balance comments: "N base difference", "N % difference",
  "N°C difference" whenever the lengths, GC % or Tm differ (by 1 or
  more); GC % and Tm carry the "!" of the single-primer criteria;
- caution and dimer flags ("!"): 3' end matches from 3 and adjacent
  homologous bases from 7 (the criteria's "reject if >="); but G/C
  runs at the 3' end and
  runs of bases are flagged from 4, one more than the criteria's
  "reject if >= 3" (a value of 3 is not flagged here, while the
  Primer Analysis screen fails it). "Runs found" and "Hairpin(s)
  found" appear when either primer is flagged; hairpins show OK or !.


## Multi Primer Analysis (recorded 2026-10-19)

Analyze -> Multi Primer asks for "primer #s (from Active Primer
List)" in four fields, Primer A to Primer D: each field takes one
digit (the list holds up to 6 primers) and moves on by itself; Enter
leaves a field empty; the screen appears after the fourth field. It
shows each primer (with its molecule position; "bp 1" for a primer
typed with Direct Entry or loaded from a file), "Primer Balance and
Cautions" (length, % GC, Tm, runs of G/C at the 3' end, runs of
bases, hairpins), and two "Primer Dimers" matrices, "3' end" and
"Any" (number of adjacent bases in the largest duplex), with the
self-dimers on the diagonal, the pairs above it, and dots below it.

`multi/multi_primer_balance.tsv` and `multi/multi_primer_dimers.tsv`:
that screen for 9 selections of 2 to 4 primers (one row per filled
field, and one row per cell of the upper triangle), among them a
field left empty (`1_2_`), the same list in another order (`2143`),
and one primer in two fields (`11__`). Recorded with
`tools/multi_primer.sh` (the primers are loaded from `.PD2` files
written by `tools/make_pd2.py`: typing four primers with Direct Entry
needs too many keys for one autotype line) and
`tools/parse_multi.py`.

Rules (`tests/unit/test_legacy_multi_primer.cpp`), the same as on the
Primer Pair Analysis screen: the balance values and their "!" are
those of "Primer pair analysis" (runs and G/C runs flagged from 4);
a diagonal cell is the primer's self-dimer values, a cell above it
the A:B values of the two primers (field order does not change
them), with flags from 3 (3' end) and 7 (any).
Not recorded: the pages behind M (molecule map: "not from molecule"
for these primers), H (hairpins of A and B, PgDn for C and D) and D
(the dimer matrices again, with the duplex of the highlighted cell),
and P (print).


### The hairpin page (recorded 2026-10-19)

`multi/multi_primer_hairpins.tsv`: the pages behind H (PgUp: A and B,
PgDn: C and D) for four primers with hairpins (from
`inspect/hairpin_design.tsv` and `inspect/inspect.tsv`), recorded with
`tools/multi_primer.sh` (`EXTRA="h , , pagedown , ,"`) and
`tools/parse_multi.py hairpins`. Each primer shows "Hairpin Energy =
N kcals" (or "No Hairpins Found") and the same drawing as the Primer
Analysis screen, one column further right (the table drops the
border column, which gives the columns of `inspect/*.tsv`). The
energies and drawings equal those of the Primer Analysis screen.

### The dimer page (recorded 2026-10-19)

`multi/multi_primer_duplexes.tsv`: the page behind D for the four
4-primer selections of `multi/multi_primer_balance.tsv`, recorded with
`tools/multi_primer.sh` (`EXTRA="d right right right down left left
down right down"`: the highlight walks A:A, A:B, A:C, A:D, B:D, B:C,
B:B, C:C, the last moves being cut by the run's time limit; then
`EXTRA="d right right right down down down"` for C:D and D:D: all 10
cells of each selection) and `tools/parse_multi.py duplexes`. The
page repeats both matrices; the highlighted cell gets "Interactions
between primers X and Y" and a drawing: the row primer 5' to 3'
(`5'-...-3'`, first base in column 11), and the column primer 3' to
5' below it (`3'-...-5'`), each base under its partner, with ':' under
each pair of the duplex.

All 40 drawings follow these rules
(`tests/unit/test_legacy_dimer_drawing.cpp`):
- the duplex is the alignment behind the "3' end" value: pairs counted
  from a 3' end, the partner extending beyond it;
- the row primer's 3' end is drawn, unless the column primer's 3'-end
  count is larger (equal counts: the row primer);
- among alignments of equal count, the one facing the partner position
  farthest from the partner's 5' end (the largest position).
When the next cell is drawn, the header changes last: the parser keeps
the first complete screen after each header.

### Hairpin drawings (from the recordings of 2026-10-17 and 2026-10-18)

The 168 hairpin drawings of `inspect/inspect.tsv`,
`inspect/hairpin_design.tsv` and `inspect/inspect_2.tsv` (columns
`hairpin_5prime`, `hairpin_pairs`, `hairpin_3prime`, leading spaces
included) follow these rules (`tests/unit/test_legacy_hairpin_drawing.cpp`):
- the 5' line shows `5'-` then the primer from its 5' end through the
  stem and the larger half of the loop; its first base is always in
  column 16;
- the marks line has `:` under each stem pair (G.T pairs included),
  then spaces, and `)` one column after the 5' line;
- the 3' line shows `3'-` then the rest of the primer, reversed, each
  base under its partner; the middle base of an odd loop ends the 5'
  line and starts the 3' line;
- nothing is drawn left of column 1: a long 3' tail loses its label
  and its 3'-most bases (one recorded case).

## 5' extension and the criteria (recorded 2026-10-18)

`inspect/inspect_extension.tsv`: four primers typed with a 5'
extension (`tools/inspect_primer.sh WORKDIR PRIMER EXTENSION`), then
Inspect, H and D; the first two columns are the extension and the
primer, then as `inspect/inspect.tsv`. The extensions make a run
across the boundary (AAAA + A...), a self-complementary site (GAATTC,
GCGGCCGC), or a hairpin.

Every criterion is computed on the whole sequence, extension
included (GC %, runs, 3' GC run, dimers, hairpins; the dimer and
hairpin drawings show the whole sequence), except the Tm, which uses
the primer without its extension (see "Primer Tm"). Values on the
Inspect screens: CCCC + GTAAAACGACGGCCAGT: GC 61 %, Tm 64 °C (as
without the extension), hairpin -2.7, 3' dimers 2, dimers 4, runs 4;
AAAA + ACGTCAGTCGATCGTCAGCT: GC 45 %, Tm 71 °C, no hairpin, 3' dimers
4, dimers 6, runs 5; GAATTC + CAGTCAGACTGACAGTCAGT: GC 46 %, Tm
59 °C, hairpin -0.5, 6, 6, runs 2; GCGGCCGC + CAGTCAGACTGACAGTCAGT:
GC 64 %, Tm 59 °C, hairpin -0.5, 6, 8, runs 2.


## Primer files `.PD2` (recorded 2026-10-18)

Written by Primer/Create -> Save Active Primer (file name typed
without extension; the program adds `.PD2`):
- `primers/M13EXT.PD2`: Direct Entry, 5' extension GAATTC, primer
  GTAAAACGACGGCCAGT (built-in criteria);
- `primers/COMP.PD2`: from molecule `pairs/TEST600.GB`, position 400,
  complement strand, 20 bases (TGTATCCAGGTAATGCAAGG).

Layout (little-endian): a 14-byte header, `1A 53 45 53` (`\x1aSES`)
then `03 01 00 00 05 00 00 00 00 00` in both files; then records of a
16-bit tag (high byte: group, low byte: field) and data, ended by tag
`FFFF` and four zero bytes.
- Group 8, text: a 32-bit length, then the text, NUL-terminated and
  padded with NUL to an even length. Fields seen: 0 = 5' extension
  (absent when there is none), 1 = primer typed with Direct Entry,
  2 = molecule name, 3 = molecule file path (`C:\TEST600.GB`), 4 =
  primer taken from the molecule.
- Group 9, numbers: 8 bytes, a zero 32-bit word then the value(s).
  Fields: 0 = length including the extension (23, 20), 1 = length
  without it (17, 20), 2 = extension length (6, 0), 3 = position on
  the molecule (1 for Direct Entry, 400), 4 = complement strand (0, 1);
  then the primer's own criteria: 0x10 = GC max, min (two int16:
  60, 50), 0x11 = Tm max, min (80, 55), 0x12 = hairpin cutoff
  (float32, 0.0), 0x13 = adjacent homologous bases, 3' dimer matches
  (7, 3), 0x14 = the two run thresholds (3, 3; their order is not
  known, both being 3).


## Uniqueness and false priming (recorded 2026-10-18)

`homology/homology.tsv`: the Uniqueness (U) and False Priming (F)
screens of primers made from molecules (Inspect, or Analyze -> Active
Primer for a complement-strand primer), with built-in settings
(exclude Tm <= 0 °C: Y; sort homologies by Tm: Y). Molecules:
`pairs/TEST600.GB` ("Primer is unique"), and `homology/TREP.GB` and
`homology/TREP2.GB`: TEST600 with planted copies of parts of primer
CTCAAAGGCTCAGACATCAG (position 51) and TGTATCCAGGTAATGCAAGG
(position 400, complement strand), on either strand, at the 3' end or
inside.

All 7 screens follow these rules:
- both strands of the molecule are searched for exact matches between
  a stretch of the primer and the molecule; the primer's own site (a
  full-length match) is not listed;
- uniqueness lists maximal matches of 7 bases or more; false priming
  lists maximal matches of 6 bases or more that include the primer's
  3' end;
- each site shows the matching bases, the molecule position facing
  the 5'-most matching base, "(C)" for a match on the complement
  strand, and the Tm of the matching bases (the primer Tm formula,
  truncated);
- sites with a Tm <= 0 are not listed (so 6-base false-priming sites
  rarely appear), and sites are sorted by decreasing Tm (sites of
  equal Tm: false-priming sites in molecule order, whatever their
  strand; uniqueness sites top strand first; see "The homology
  settings" below).
Not probed yet: the "exclude Tm <= 0" and "sort by Tm" settings set
to N, and molecules with more sites than one screen holds.


### The homology settings (recorded 2026-10-19)

`homology/homology_settings.tsv`: Uniqueness and False Priming screens
(in screen order: `rank`) with the System Defaults "False priming and
uniqueness search: exclude Tm <= 0" and "Sort homologies by Tm?" set
to Y or N. The settings come from files loaded as `PD2.CFG`:
`config/PD2_EXCLUDE_N.CFG`, `config/PD2_SORT_N.CFG` and
`config/PD2_EXCLUDE_SORT_N.CFG` are not recordings, but inputs written
by `tools/make_pd2_cfg.py` (built-in defaults and the changed flags;
the tool reproduces `PD2_LENGTH25.CFG`, `PD2_WEIGHTS12345.CFG` and
`PD2_WEIGHTS35142.CFG` byte for byte). Recorded with
`tools/homology_screen.sh` (`CFG=...`) and `tools/parse_homology.py`.
Molecules: `homology/TREP2.GB`, and three written for these probes
with `tools/plant_sites.py` (the planted copies are in their
DEFINITION line): `homology/THOM.GB` (copies of 6 to 9 bases of the 3'
end of the primer at 51, GGGGTCGAGCAAGACCGTAT, on both strands),
`homology/THOM2.GB` and `homology/THOM3.GB` (two copies of the 3' end
of the primer at 51, AGCAATTTGATCATGCAGAT, on both strands; in
`THOM3.GB` both copies give sites of 9 bases and Tm 10).

All 17 screens follow these rules:
- exclude Tm <= 0 set to N: sites whose Tm is not positive are listed
  too, with Tm 0 (6- and 7-base false-priming sites, for instance);
- sort by Tm set to N: false-priming sites are listed by molecule
  position (the position shown), whatever their strand; uniqueness
  sites are listed by strand (top strand first), then by position;
- sort by Tm set to Y: that order, sorted by decreasing Tm, sites of
  equal Tm keeping it. So false-priming sites of equal Tm are in
  molecule order, whatever their strand (`THOM3.GB`: 208 (C) before
  400), as the section above states.

## Match Active Primer (recorded 2026-10-18)

`search/match_active.tsv`: the Find Matching Primer list for primer A
at position 51 of `pairs/TEST600.GB` (CTCAAAGGCTCAGACATCAG, GC 50 %,
Tm 66 °C), with the form's defaults: target 400, search from target
- 50 to target + 50, complement strand, length 20, no extension, GC
within 5 % and Tm within 10 °C of primer A. Keys: as for a molecule
primer, then s, down down enter (Match Active Primer), 400 enter, and
enter for the 7 other fields; the list was scrolled with up and down
to see all 11 entries.

The list holds exactly the primers B (20 bases, complement strand,
position = 5' end, from 350 to 450) that:
- meet every single-primer criterion (see "Primer criteria");
- have a GC % within 5 and a Tm within 10 °C of primer A (displayed
  values);
- have fewer A:B 3' end matches than the criterion's "reject if >="
  (3); 371, with 5 adjacent homologous bases A:B, is listed, so the
  A:B adjacent bases (all below 7 here) were not tested by this
  recording.
Entries are sorted by position; "bps from target" is the distance
between the position and the target. Every primer of the window that
passes these tests is listed, and no other.


## Find Pair (recorded 2026-10-18)

`search/find_pair.tsv`: the Find Primer Pair list for
`pairs/TEST600.GB` with targets 150 (primer A) and 400 (primer B) and
the form's other defaults: primer A searched from target - 100 to the
target, primer B (complement strand) from the target to target + 100,
length 20, no extensions, GC within 5 % and Tm within 10 °C, all five
weights 1. Keys: space, down down down enter (Find Pair), enter, down
down enter (TEST600), enter (Bp position), 150, enter x3, 400, enter
x5. The program reported "Primer A: 3 possible primers, Primer B: 7
possible primers, Pairs: 10 possible", then the list (rank = order on
screen).

Rules that reproduce the list:
- possible primers: those of the search windows that meet every
  criterion (A: 71, 72, 73; B: 413, 414, 422, 424, 425, 453, 457);
- possible pairs: GC % within 5, Tm within 10 °C, and fewer A:B 3'
  end matches than the criterion (3), as in Match Active Primer;
- "Max Dimers 3'" and "Any": the largest of the A:A, B:B and A:B
  values; "FalsePriming TmA/TmB": the highest false-priming site Tm
  of each primer on the molecule, "-" if none;
- order: with all weights 1, the eight pairs without false-priming
  site are sorted by GC difference + Tm difference + max 3' dimers +
  max adjacent bases (ties: primer B's position); the two pairs whose
  primer A has a false-priming site (Tm 21) come last.

`search/find_pair_weights.tsv`: the same search with other weights
(set with F4 in the form; each weight field takes one digit and moves
to the next field by itself, so no Enter is typed between them). The
same 10 pairs are found; their order:
- weights 0 0 0 0 1 (false priming only): the eight pairs without a
  false-priming site first, in generation order (A, then B, by
  position), then the two with one;
- weights 0 5 0 0 1: by 5 x Tm difference, ties in generation order;
  the pairs with a false-priming site come after pairs of score 45,
  their own score being 30 or 45 plus the false-priming term.

So the order is a stable sort, by increasing score, of the pairs in
generation order, with score = sum of weight x term, the terms being
GC difference, Tm difference, max 3' dimers and max adjacent bases;
the false-priming term (weight 1) is above 15 for a site of Tm 21.
Taking the site's Tm itself fits all three lists, but is not proven.

### Limits on primers and pairs (probed 2026-10-18, partial)

`config/PD2_RELAXED.CFG` is not a recording: it is a settings file
written for probes, with the layout of "Settings file `pd2.cfg`"
(built-in values, except GC 0-100 %, Tm 0-99 °C, hairpin cutoff -99
kcal, dimer and run thresholds 99), so that nearly every primer meets
the criteria. Copied as `PD2.CFG` next to PRIMER.EXE, it is loaded at
start (a file with another name is ignored). Find Pair on
`enzymes/TLONG.GB` (1200 bases, random), targets 400 and 800:
- windows of 41 positions (target - 40, target + 40): "Primer A: 41
  possible primers", "Primer B: 41", "Pairs: 626 possible", then the
  list;
- windows of 61 positions: 61 and 61 possible primers; the pair count
  climbs to "Pairs: 1000 possible" and stops there, while the search
  goes on ("Press key to interrupt search") for more than 20 minutes
  of emulated time, without another message. Whether a warning
  follows ("More than 1000 pairs", proposal §1.2) was not observed.
- the 100-site message ("More than 100 sites found. Suggest you
  redefine matching limits.") was not reached: the run with windows
  of 151 positions lost keys. To probe again with long time limits,
  one emulator at a time.
Keys typed during the search may interrupt it: type nothing after the
last Enter of the form until the counts are final.

### Limits on primers and pairs (recorded 2026-10-19)

`search/find_pair_limits.tsv`: Find Pair on `enzymes/TLONG.GB` with
chosen search windows, with `config/PD2_RELAXED.CFG` as `PD2.CFG`
(nearly every primer meets the criteria) or the built-in criteria.
Recorded with `tools/find_pair_limits.sh` (one emulator at a time,
`PAUSE=4`; runs that lost keys were discarded). Offsets are those of
the form ("Start search at: Target -", "Extend search to Target +");
a window holds offset + 1 positions. "-": not shown (the message
comes while the primers are counted). Outcomes:
- `pairs_found`: the counts are final and "<Press key to display
  possible pairs>" follows;
- `too_many_sites`: while primer A, then primer B, are counted, a box
  shows "More than 100 sites found.  Suggest you redefine matching
  limits." and "Press any key to continue..."; a key closes it and
  the constraints form is shown again, with no search;
- `pairs_capped`: the pair count climbs to "Pairs:  1000 possible"
  and stops; "<Press key to interrupt search>" stays on screen for
  more than 40 minutes of emulated time, and keys (space, Esc, y,
  Enter, sent 1 minute after the cap) change nothing: the program
  hangs. A control run (41x41 windows, same delayed key) shows the
  list, so the keys do reach the program;
- `no_matching_primers`: a box "No Matching Primers Found" and "Press
  any key to continue...", without the search screen.

Rules:
- the limit is on possible primers (those that meet every
  criterion), not on window positions: 111 positions with 9 possible
  primers A (built-in criteria) are searched; it applies to primer A
  and to primer B;
- the message appears at 100 possible primers, despite its wording:
  99 possible primers A (relaxed, offset 98) are searched, 100
  positions (offset 99, every primer possible with these settings)
  give the message;
- the pair count stops at 1000. The program holds the string "More
  than 1000 possible pairs found!" but never showed it.
- the hang does not happen every time: on 2026-10-20, Find Pair on
  `pairs/TEST600.GB` (targets 150 and 400, offsets 100, with
  `config/PD2_CRITERIA.CFG`: 37 possible primers A, 45 B) also stopped
  at "Pairs:  1000 possible", then showed its list and printed it
  (`reports/find_pair/TEST600_150_400_CRITERIA.PRN`, see "The Find Pair
  report"); what makes the recorded runs above hang is not known.
Not known: what makes "No Matching Primers Found" appear (no
possible primer A? no pair?); in both recorded cases no pair exists.
Version 1 (decided 2026-10-19): the same 100-primer limit; at 1000
pairs, the search stops, keeps the first 1000 pairs in generation
order, ranks them, and warns "More than 1000 possible pairs found!"
instead of hanging.

### No Matching Primers Found (recorded 2026-10-19)

`search/find_pair_no_match.tsv`: Find Pair on `enzymes/TLONG.GB` with
built-in settings and the form's windows (100 and 100), recorded with
`tools/find_pair_limits.sh`, for targets 120/300 (primers A but no
possible primer B, as primerpicker computes them) and 500/1100
(possible primers on both sides, but no pair within the GC and Tm
tolerances). With the two `no_matching_primers` rows of
`search/find_pair_limits.tsv` (no possible primer A), all three cases
give the same box, "No Matching Primers Found" and "Press any key to
continue...", without the search screen or any count. Rule: the
message replaces the list whenever no pair is found. (The CLI keeps
its count lines, which an existing test pins, and prints the message
instead of the empty list.)
Navigation: a target field takes at least 6 digits (it does not move
on by itself after 4).

### The false-priming term of the score (recorded 2026-10-18)

`search/find_pair_false_priming.tsv`: Find Primer Pair lists of
`search/TFP.GB` and `search/TFP2.GB` (targets 150 and 400, other
settings built in), with several weight sets. Both molecules are
random sequences with planted copies of primer 3' ends
(`tools/plant_sites.py`; the copies are listed in the DEFINITION
line), which give false-priming sites to primers A (site Tm 15 to
39) and to one primer B (site Tm 6); in `TFP2.GB`, seven primers A
have a second, weaker site. 30 pairs per list. Recorded with
`tools/find_pair.sh` (`DOWNS=32`, `PAUSE=2.5`) and
`tools/parse_find_pair.py`. Runs that lost keys were discarded:
the parser checks the first four weights on the weights form (the
form closes before showing the fifth) and that the list was
reached.

All 12 recorded lists (these 9, and the 3 lists of `TEST600.GB`)
follow this score, by increasing score and in generation order for
equal scores:

    score = wGC x |GC% A - GC% B| + wTm x |Tm A - Tm B|
          + w3' x max 3' dimers + wAny x max dimers (any)
          + wFP x (false-priming Tm of A + false-priming Tm of B)

- the false-priming Tm of a primer is that of its highest site (the
  TmA/TmB columns), 0 when it has none; the sum over all its sites
  does not fit the `TFP2.GB` lists;
- the term is the sum of both primers' Tm: the largest, or primer
  A's alone, contradict several lists;
- the weight multiplies it (`0 5 0 5 2`: the order differs from
  `0 5 0 5 1` predictions);
- a pair whose only site has Tm 6 is placed exactly as a score of 6.

Not probed: lists of more than 1000 pairs or more than 100 sites
(the program warns), and other targets or windows.


## Mutagenize (recorded 2026-10-19)

`mutagenize/mutagenize.tsv`: the Mutagenize Primer screen (Modify ->
Mutagenize) of 20-base primers made from `pairs/TEST600.GB` at
positions 5, 51, 52, 53, 200 and 575, two of them after edits,
recorded with `tools/mutagenize.sh` (KEYS typed on the screen) and
`tools/parse_mutagenize.py` (the lines inside the screen's box, the
border column removed). Editing: the primer is typed over (overwrite
mode); a typed base replaces the one under the edit highlight, which
does not move; the arrow keys move it.

All 6 screens follow these rules (`tests/unit/test_legacy_mutagenize.cpp`):
- the molecule line shows the molecule from 20 bases before the
  primer to 39 bases after its 5' end, cut at the molecule's ends; the
  primer's first base is in column 27, under its molecule base;
- three rows above (labelled O, R, I: "ORI"ginal) translate the
  molecule with the standard genetic code (stop codons as `-`): each
  amino acid under the third base of its codon, codons entirely
  inside the molecule window, starting at positions of 0, 2 and 1
  modulo 3 (O, R, I; 1-based positions);
- three rows below the primer (M, U, T: "MUT"ated) translate the
  edited primer: codons entirely inside the primer, starting at
  positions of 1, 2 and 0 modulo 3;
- the ORIGINAL and MUTATED lines give %GC and Tm (with `!` outside the
  criteria), then Hairpin, Dimers (both dimer criteria) and Runs (both
  run criteria) as OK or `!`;
- the line under the molecule always reads `5'->3'` and `--` (in
  columns 32-33, at every position recorded).
`mutagenize/codon_table.tsv`: the codon table behind F2 (a scrolled
list, opened on the amino acid under the edit highlight; `home` goes
to the top), recorded with `KEYS="f2 , home , pagedown ..."`: the 20
amino acids by three-letter name, with their one-letter code and
codons (alternatives for the third base in parentheses), no line for
stop codons. It is the standard genetic code
(`tests/unit/test_codon_table.cpp` expands every line).
Not recorded yet: F1 (the enzyme list), S (save), I (inspect), a
complement-strand primer, a primer with a 5' extension.

## Reports (recorded 2026-10-19 and 2026-10-20)

`reports/*.PRN`: the reports PRIMER.exe prints, written with P (or the
Print menu of the Uniqueness and False Priming screens) -> Print
Control -> Print to Disk, file name REPORT (the menu also offers
Abandon Printing and Start Printing); `reports/reports.tsv` lists the
screen and inputs of each file. Recorded with `tools/print_report.sh`
(and `tools/multi_primer.sh` for the Multi Primer report). The files
are kept as written: DOS line ends, code page 437 (the degree sign is
byte F8); PAIR.PRN has CR and a form feed (bytes 0D 0C) before its page
2, the other reports none.

What they hold (One Primer, Primer Pair, Multiple Primers and Find Pair
are now reproduced by primerpicker, see the sections below; Uniqueness
and False Priming are not yet):
- every report starts with a title line ("Primer Design Report:  One
  Primer", "...:  Multiple Primers", "...:  Find Pair", "Uniqueness of
  Primer Report", "False Priming Report") and "Report Date: MM/DD/YY"
  (the computer's date);
- INSPECT.PRN (Primer Analysis): the primer, its length, 5' extension
  and composition, the criteria table (setting, actual value, meets
  criteria), then "LARGEST DUPLEXES REPORTED IN PRIMER DIMERS": the
  self-dimer with 3' end involvement, and the other (longest) self-
  dimer, drawn as on the dimer page but with the first base in
  column 13;
- PAIR.PRN (Primer Pair Analysis) and MULTI.PRN (Multi Primer): the
  primers, the balance and cautions, the dimer matrices; the Multi
  report prints % GC without the '!' of the screen (45 for primer B);
- FINDPAIR.PRN: the source molecule and its size, the search
  parameters, the criteria, then the best pairs (P first asks "N
  Pairs found.  Print best: [N]") in the order of the list;
- UNIQUE.PRN and FALSEPRI.PRN: the primer and its position, the
  source molecule, the search criteria (with "Sites with melting
  temperatures <= 0 C are not reported"), a text map of the molecule
  (`P>` at the primer, `>` and `<` for sites on the top and
  complement strands), and the sites as on the screen.

### The One Primer report (recorded 2026-10-20)

`reports/inspect/*.PRN`: 46 more reports of the Primer Analysis screen,
written with `tools/print_report.sh inspect PRIMER [EXTENSION]` (with
`CFG=config/PD2_....CFG` for a settings file); `reports/inspect_reports.tsv`
lists every One Primer report (INSPECT.PRN included): file, primer,
5' extension, settings file, and whether both duplexes have pairs.
Primers: those of `inspect/` with ties between duplexes, without
complementary bases, with hairpins, of 8 to 46 bases; four
extensions (1 to 12 bases); `config/PD2_RELAXED.CFG`,
`config/PD2_REPORT_NARROW.CFG` (GC 40-65 %, Tm 50-70 °C, hairpin cutoff
-3.5 kcal, dimers 4 and 6, runs 4, 3' G/C runs 2) and
`config/PD2_REPORT_WIDE.CFG` (GC 5-100 %, Tm 100-120 °C, cutoff -12.5
kcal, dimers 10 and 12, runs 9 and 9), written with
`tools/make_pd2_cfg.py --criteria` (inputs, not recordings).

Every report with pairs is reproduced byte for byte by
`legacy_inspect_report` (`tests/unit/test_legacy_inspect_report.cpp`,
`tests/cli/test_report.sh`):
- the sequence is printed `5'-EXTENSION-PRIMER-3'` (the extension, a
  dash, the primer); Length counts both; "5' Extension:" shows the
  extension or `--none--`;
- the criteria rows: the setting starts in column 17, the actual value
  ends in column 56 (a long setting, such as `-99.0 kcal`, pushes
  nothing: the value keeps its right edge), then 11 spaces and YES or
  NO; the hairpin row shows `--`, as on the screen; no trailing spaces;
- the two duplexes, see "Duplexes of the reports" below.

**Duplexes without pairs.** When a duplex has no pair (the 3' end of a
primer pairs with nothing, or the primer has no complementary bases at
all: AAAAAAAAAAAA, ACACACACACAC, TGTGTGTGTGTGTGTG, ACGCCGTTTTCGGTGT,
ACTCCCTTTTGGGGGT), PRIMER.exe writes lines of tens of thousands of
spaces: for the 3'-end duplex, a marks line of 51121 spaces then a line
`-5'`; for the other duplex, a marks line of 6990 (12 bases) or 7006
(16 bases) spaces and a bottom line of about 29000 spaces before
`3'-...-5'`. These reports (`duplexes` = `no-pairs`) are recorded, and
not reproduced: decided with the developer on 2026-10-20, primerpicker
draws the partner right under the primer (its `3'-` under the `5'-`)
below an empty marks line (`tests/unit/test_report_no_pairs.cpp`).

### The Primer Pair report (recorded 2026-10-20)

`reports/pair/*.PRN`: 43 more reports of the Primer Pair Analysis
screen, for the pairs of `pairs/pair_screens.tsv` and
`pairs/pair_screens_2.tsv` and a few others, written with
`tools/print_report.sh pair` (Direct Entry) or `pair_files` (the
primers loaded as `.PD2` files, for pairs too long to type in one
autotype line); `reports/pair_reports.tsv` lists every Primer Pair
report (PAIR.PRN included). Every report with pairs is reproduced byte
for byte by `legacy_pair_report` (`tests/unit/test_legacy_pair_report.cpp`,
`tests/cli/test_report.sh`):
- two pages: page 1 ends with an empty line, then CR and a form feed
  (bytes 0D 0C) start the title of page 2 ("-- Page 2");
- each primer: "% GC:" and "Tm C:" values in 3 columns (`% GC:100;`),
  the composition without padding, the criteria in short form, "Meets
  criteria selected: YES/NO" (as the Primer Analysis screen);
- Primer Balance: no '!' flag; the comment ("N base difference", "N %
  difference", "N C difference") only when the difference is not zero,
  6 spaces after the value of B;
- Primer Cautions: values followed by their '!' flag or a space (hence
  trailing spaces), as the screen; Hairpins shows OK or '!' with no
  flag column, then "Hairpin(s) found" (5 spaces after B) when either
  primer has one; no "Runs found" comment, unlike the screen;
- Primer Dimers: A:A, B:B, A:B with flags, then "(see page 2" and "for
  duplexes)";
- page 2: the 3'-end duplexes, then the other duplexes, of A:A, B:B
  and A:B (see "Duplexes of the reports"; for A:B, A is the first
  primer).
AAAAAAAAAAAA with CCCCCCCCCCCC (no pair at all) writes the long lines
of spaces described above (`duplexes` = `no-pairs`; not reproduced).

### The Multiple Primers report (recorded 2026-10-20)

`reports/multi/*.PRN`: 6 more reports of the Multi Primer Analysis
screen, for the selections of `multi/multi_primer_balance.tsv`
(`tools/multi_primer.sh` with `EXTRA="p , , down down enter , , r e p o
r t enter , , , , , ,"`; the run keeps the report as `c/REPORT`);
`reports/multi_reports.tsv` lists every Multiple Primers report (fields
A to D, '-' for an empty field). All are reproduced byte for byte by
`legacy_multi_primer_report` (`tests/unit/test_legacy_multi_primer_report.cpp`,
`tests/cli/test_report.sh`):
- one section per field, A to D; an empty field's section holds
  `     none` instead of its sequence, source molecule and length;
  "Source molecule:  none" for primers that do not come from a molecule
  (all recordings so far);
- Primer Balance and Cautions: the label, then four columns of 8
  characters, values right-aligned, no '!' flag (the hairpin '!'
  replaces OK); an empty field leaves its column blank, trailing spaces
  included;
- Primer Dimers: the "3' end" and "Any" matrices side by side, values
  with their '!' flag, dots below the diagonal; an empty field leaves
  its row (no dots) and its column blank; these lines have no trailing
  spaces;
- no duplex drawings.
The report of selection `12__` was not written (keys lost).

### The Find Pair report (recorded 2026-10-20)

`reports/find_pair/*.PRN`: 4 more reports of the Find Primer Pair
list (`tools/print_report.sh find_pair`, with `CFG` for a settings
file: the file is listed before the molecule, which takes one more
down): `search/TFP.GB` and `search/TFP2.GB` (30 pairs each,
false-priming Tm values),
`pairs/TEST600.GB` with `config/PD2_CRITERIA.CFG` (1000 pairs: the pair
limit; this search did not hang, and its best 20 pairs are those of the
first 1000 pairs, ranked, as primerpicker keeps them) and with
`config/PD2_WEIGHTS35142.CFG` (10 pairs); `reports/find_pair_reports.tsv` lists
every Find Pair report (FINDPAIR.PRN included). All are reproduced byte
for byte by `legacy_find_pair_report`
(`tests/unit/test_legacy_find_pair_report.cpp`, `tests/cli/test_report.sh`):
- P first asks "N Pairs found.  Print best: [M]", M being N up to 20
  (30 pairs: `[20  ]`); the recordings keep M (Enter): the report
  holds the best M pairs, in the order of the list;
- "Source molecule:" then the molecule's name, then its size
  right-aligned to end in column 56, then "bps" (see "The Uniqueness
  and False Priming reports" for sizes of 4 and 5 digits);
- the search: targets in 6 columns, the offsets, the primer length (in
  4 columns after "Length:"), "5' Extension: - none -", the GC and Tm
  tolerances; the weights are not printed;
- the criteria in short form, as the Primer Pair report ("Hairpin
  cutoff:" then the cutoff in 4 columns: `cutoff:-1.5`);
- the header of the list has the degree sign of code page 437 (byte
  F8); the rows are those of the screen, a false-priming Tm of B taking
  2 columns (`- /  6`), '-' when there is no site.
Not recorded: reports of more than 20 pairs (whether PRIMER.exe breaks
pages), a typed M. `config/PD2_LENGTH25.CFG`
lost keys twice.

### The Uniqueness and False Priming reports (recorded 2026-10-20)

`reports/homology/*.PRN`: 21 more Uniqueness and False Priming reports
(`tools/print_report.sh homology MOLECULE POSITION u|f`, with `CFG`
for a settings file), on molecules written for these probes with
`tools/make_map_molecule.py` (`homology/map/`, inputs, not recordings):
random sequences of 100 to 50,000 bases with copies of a 12-base
stretch of the primer at 51 (55:66, or 59:70 with the 3' end for
`TMAP660F.GB`) at chosen positions, on the top strand (t) or as
reverse complements (c), as listed on their DEFINITION line. Primers at
positions 1 to 900; settings files `PD2_LENGTH25.CFG`,
`PD2_EXCLUDE_N.CFG`, `PD2_SORT_N.CFG`. `reports/homology_reports.tsv`
lists every report (UNIQUE.PRN and FALSEPRI.PRN included).

The molecule map (all 23 reports, 122 sites):
- one character stands for s = L / 73 + 1 bases (integer division; L
  is the molecule's length): 100 bases: 2; 600: 9; 660: 10; 730: 11
  (not 10); 1000: 14; 6600: 91; 50,000: 685;
- the line of `=` has L / s characters (50, 66, 71, 72, ...);
- the line above it holds `P>`, starting in column max(0, P / s - 1),
  P being the primer's position;
- below it, a line for the top strand (`>`) and one for the complement
  strand (`<`): each site is marked in column 3 + (Bp Pos - 1) / s,
  Bp Pos being the position the site's row shows; sites in the same
  column give one mark (`<<<` is three columns); a line without marks
  is a single space.
The rest of the report:
- title, "Primer Sequence:", "Bp pos:" (6 columns), "Source
  molecule:" (the size right-aligned, ending in column 56: `600`,
  `1000`, `50000`), the criteria ("Regions of homology (>= 7 bases)"
  or "3' end homologies (>= 6 bases)"), then "Sites with melting
  temperatures <= 0 C are not reported", which is an empty line when
  that setting is N;
- one row per site, in the order of the screen: the matching bases
  right-aligned after dots in a field of 19 characters (whatever the
  primer length: 20 and 25 recorded), then `....3'` for a uniqueness
  site or 6 spaces for a false-priming site, the position ending in
  column 52, ` (C)` or 4 spaces, the Tm ending in column 66, and 3
  spaces;
- when there is no site, PRIMER.exe shows a box ("No false priming
  found") and prints nothing.
All are reproduced byte for byte by `legacy_homology_report`
(`tests/unit/test_legacy_homology_report.cpp`), except one:

**Primers with an internal repeat.** `homology/map/TOWN.GB` and
`homology/map/TOWN2.GB` (`tools/plant_sites.py`, 1000 bases each) hold
20-base primers with a direct repeat of r bases at primer offsets a
and a + d (0-based), an inverted repeat, and a control without repeat
(TOWN at 801); `TMAP1000.GB` at 900 is one more case. PRIMER.exe's
Uniqueness lists differ from a search of maximal matches for every
primer with a direct repeat of 7 bases or more (`note` =
`own-site-overlap`; not reproduced, primerpicker lists the maximal
matches):

| primer | r | a, d | listed: primer offsets, at | maximal matches not listed |
|---|---|---|---|---|
| TMAP1000 at 900, AAATTCAAAATTCATGTGAT | 7 | 0, 7 | 1-12, at 901 | |
| TOWN at 201, ACGTTGCTAACGTTGCGAGT | 7 | 0, 9 | 1-14, at 202 | AACGTTGC at 200 |
| TOWN at 351, GATCCAGTACGATCCAGTTT | 8 | 0, 10 | 2-15, at 353 | CGATCCAGT at 350 |
| TOWN at 501, TGCAGGTCACTGCAGGTCAA | 9 | 0, 10 | 3-15, at 504 | TGCAGGTCA at 501 and 511 (Tm 17) |
| TOWN at 651, CCTAGTAGATTGTCCTAGTA | 7 | 0, 13 | 1-18, at 652 | |
| TOWN2 at 351, GTCAGGACTTGTCAGGACTT | 10 | 0, 10 | 4-15, at 355 | TGTCAGGACTT at 350, GTCAGGACTT at 361 |
| TOWN2 at 651, GATTCCGATTCAGATTCCGA | 8 | 0, 12 | 2-17, at 653 | GATTCCGAT at 663, GATTCCGA at 651, ATTCCGAT at 703 (C) |
| TOWN2 at 801, TCAGCTTAGCAGCTTAGCTC | 9 | 1, 8 | 4-14, at 805 | CAGCTTAGC at 802 and 810 |
| TOWN2 at 201, (ACGTC) x 4 | 15 | 0, 5 | none: "Primer is unique" | four matches of 10 and 15 bases |

- the alignments of one copy of the repeat on the other, inside the
  primer's own site, are not listed (even with a positive Tm);
- a stretch of the primer's own alignment is listed instead: from
  offset a + r - 6 to offset a + d + 5 (the last 6-base window of the
  first copy to the first 6-base window of the second), when it has 7
  bases or more; all 8 cases above (the tandem repeat gives 2 bases:
  nothing). Predicted before recording for the four TOWN2 cases;
- a repeat of 6 bases changes nothing (TOWN2 at 51,
  GACCTATCGACCTACTGAGT, r = 6, d = 8: "Primer is unique", although the
  rule above would give 14 bases); an inverted repeat neither (TOWN2 at
  501: "Primer is unique", as primerpicker);
- a site outside the own site can vanish too: ATTCCGAT at 703 (C), Tm
  9, for TOWN2 at 651; its 6-base windows ATTCCG and TTCCGA occur
  twice in the primer;
- False Priming is not affected: TOWN2 at 651 lists GATTCCGA at 651
  (the second copy, at the 3' end, on the first), as primerpicker does
  (`reports/homology/TOWN2_651_f.PRN`, reproduced byte for byte); TOWN
  at 651 and TOWN2 at 801: "No false priming found", as primerpicker.
This suggests a search by 6-base windows whose own-site and
duplicate-window handling goes wrong for windows that occur twice in
the primer; no model tried so far reproduces every case.
Version 1 (decided 2026-10-20): primerpicker keeps its search of
maximal matches, which lists every real match; these recordings
document a known deviation from PRIMER.exe.

Not recorded: a complement-strand primer (the tool makes top-strand
primers only); a primer within 20 bases of the molecule's end (two
runs froze); "Primer is unique" (the uniqueness box, presumably
without report either).

### Duplexes of the reports (from the recordings of 2026-10-20)

The One Primer report ends with "LARGEST DUPLEXES REPORTED IN PRIMER
DIMERS": the primer (with its extension) with itself, drawn as on the dimer page (the first
primer 5' to 3', first base in column 13; the second 3' to 5' below
it; ':' under the pairs of the duplex):
- "Primer-dimers with 3'end involvement": the duplex of the dimer page
  (see "The dimer page"); its length is the "3' Dimers" value;
- "Other Primer-Dimers": the longest run of consecutive pairs in any
  alignment (its length is the "Dimers" value); among equally long
  runs, the first one in the order of the alignments, from both 5'
  ends facing each other on (second primer shifted to the right), then
  from the first primer's 5' end. All 42 reports with pairs follow
  this rule; most of them have tied runs to choose from.
- the second primer's line: `3'-` is left out when it does not fit
  left of the 3' base (the bases keep their columns), and bases left
  of column 0 are cut (`AGTTTCTTTTGAGGCT-5'` with one leading space,
  `ATCCGGCGAAGTTAGCT-5'` from column 0).

## Circular molecules (recorded 2026-10-20)

`circular/`: two pairs of molecules written with `tools/plant_sites.py`
(copies that run past the last base continue at the first), each
identical but for the LOCUS line: `TCIRC.GB` (circular) and `TLIN.GB`
(linear), 600 bases, with a primer meeting the built-in criteria across
the origin (CTCAATGACGTGCGGTTAGT, bases 591 to 600 and 1 to 10) and a
complement-strand primer B ending at 50; `TCIRC2.GB` and `TLIN2.GB`,
with the primer CTCAAAGGCTCAGACATCAG at 300 and copies of its bases 4
to 15 at 450 and across the origin (595 to 600 and 1 to 6). Recorded
with `tools/print_report.sh` (`find_pair`, `homology`).

Apart from the restriction site search (see "True and almost sites"),
PRIMER.exe treats a circular molecule as linear:
- Find Pair, targets 15 and 560 (offsets 100): the same list on both
  molecules; the search windows stop at the molecule's ends (primer A
  from 1, primer B up to 600), so neither the primer across the origin
  nor the primer B at 50 is found; the report prints the offsets cut
  at the ends: "Start: Target - 14", "End: Target + 40"
  (`reports/find_pair/TCIRC_15_560.PRN`, `TLIN_15_560.PRN`);
- Uniqueness of the primer at 300: the copy across the origin is not
  found on the circular molecule either, only the copy at 450
  (`reports/homology/TCIRC2_300_u.PRN`, `TLIN2_300_u.PRN`, reproduced);
- a primer made from the molecule at 591 (20 bases in 600) is refused
  with "Invalid Position", on both molecules.
primerpicker v1 does the same.

**Primers B ending with a run of three.** In these Find Pair searches,
PRIMER.exe found 9 possible primers B and 25 pairs, where a Find Pair
that applies the Primer Analysis criteria finds 11 and 33: the two left
out, at 569 (AGGCCAGCGATCACGTCAAA) and 598 (TCATTGAGCTCCTCGCGAAA), end
with AAA. See "Find Pair and runs at the 3' end" below; both reports
are reproduced.

## Find Pair and runs at the 3' end (recorded 2026-10-20)

`search/TRUN.GB` (`tools/plant_sites.py`, 2000 bases) holds, in five
regions 400 bases apart, a primer A (top strand, at 400 i + 50) and a
primer B (complement strand, ending at 400 i + 219) that meet every
Primer Analysis criterion and end with AA (region 0), AAA, AAAA, TTT and
TTTT (regions 1 to 4); in regions 1, 3 and 4 primers A and B have the
same sequence. Find Pair with targets 400 i + 80 and 400 i + 180
(offsets 100), with built-in settings and with `config/PD2_RUNS4.CFG`
(runs criterion 4, otherwise built in; written with
`tools/make_pd2_cfg.py --criteria`, an input); recorded with
`tools/print_report.sh find_pair`. `search/find_pair_3prime_runs.tsv`:
the final "Primer A/B: N possible primers" and "Pairs: N possible"
counts of the 9 runs that reached them (the counts climb during the
search: the last value is kept); `reports/find_pair/TRUN_*.PRN`: the 8
reports.

Rule (all 9 primer A counts, 7 primer B and pair counts, and
`circular/TLIN.GB` and `TCIRC.GB`): Find Pair rejects a primer whose
run at the 3' end reaches the runs criterion ("reject >= 3 base runs":
AAA is rejected, AA is not; with 4: AAAA is rejected, AAA is not),
although the Primer Analysis screen does not count that run ("Primer
criteria"). A fixed limit of 3 contradicts the counts with criterion 4.
Predicted before recording from the TLIN case, for both settings.
`tests/unit/test_legacy_find_pair_3prime_runs.cpp` checks the counts;
5 TRUN reports and both circular reports are reproduced byte for byte.

**Primers with a copy elsewhere.** Where primers A and B have the same
sequence (regions 1, 3 and 4), some primers have a second full-length
(or 19-base) match in the molecule, at the other primer's site: there,
PRIMER.exe's lists differ (two runs: the primer B count, the pair count;
three reports: the order of the pairs). primerpicker's homology search
leaves out every full-length match as the primer's own site, which
these recordings suggest is not what PRIMER.exe does (see "Copies of
the primer elsewhere"; `note` = `full-length-copy` in
`find_pair_3prime_runs.tsv`; reports in
`reports/find_pair_reports_pending.tsv`).

**Match Active Primer** applies the same rule: the Find Matching Primer
lists of primer A at 50 (CTCAATGACGTGCGGTCGAA) for targets 1020 and
1820 of `search/TRUN.GB` (`tools/print_report.sh match`, form
defaults), printed to disk (`reports/match/TRUN_50_*.PRN`), leave out
the primers B at 1019 (AAAA), 1020 (AAA), 1819 (AAAA), 1820 and 1860
(AAA), as predicted; `search/match_3prime_runs.tsv` lists the positions
(`tests/unit/test_legacy_match_3prime_runs.cpp`). The format of the
Match Active Primer report is recorded, not reproduced yet.

## Restriction enzymes (recorded 2026-10-18)

### Built-in enzyme table

`enzymes/enzyme_table.tsv`: the enzyme table compiled into PRIMER.exe,
extracted with `tools/extract_enzymes.py` (the executable is packed
with Microsoft EXEPACK; the tool unpacks it). It is recorded here to
compare search results with PRIMER.exe only: PrimerPicker ships
enzymes from REBASE, not this table (proposal §3.6).

- 359 enzymes in 28-byte records, in alphabetical order: name,
  recognition sequence (135 with IUPAC letters or N), cut position,
  overhang, length, and two flag bytes.
- Cut and overhang were checked on a few well-known enzymes (AatII
  GACGT^C: cut 5, overhang -4; EcoRI G^AATTC: 1, +4; AluI AG^CT: 2,
  0; BglI GCCNNNN^NGGC: 7, -3). They are not shown on the enzyme
  screens.
- Flag 0x80 marks 89 enzymes: the first ten (AatII, AflII, AluI, AlwI,
  ApaI, ApaLI, AscI, AseI, Asp718, AsuII) are exactly the first
  screen of the Find Enzymes menu, which shows name and recognition
  sequence. The other 270 are presumably the "ALL" list from which
  enzymes can be added (not probed yet).

### Find Enzymes

Molecule-Ops -> Find Enzymes asks for a molecule (Load File), then
offers "Almost Sites" (1 base mismatches to the recognition sequence)
and "True Sites" (exact matches), then the enzyme menu (arrows or
initial letter; PgDn does nothing). Keys are lost while the molecule
map is drawn: type them 3 s apart.

True Sites for AluI (AGCT) on `pairs/TEST600.GB`: "Number of Sites =
1", position 322, the first base of the only AGCT of the top strand.

### True and almost sites (recorded 2026-10-18)

`enzymes/find_enzymes.tsv`: True Sites and Almost Sites screens on
`enzymes/TENZ.GB`, a random linear molecule of 600 bases (seed 4015,
no site of the probed enzymes) with planted sites: AluI AGCT at 1,
150 and 597 (the first and last bases); AlwI GGATC at 200 and its
reverse complement GATCC at 300; FokI GGATG at 400 and 592 and CATCC
at 6 and 460 (FokI cuts 9/13 bases outside its site: at 6 and 592
the cut is outside the molecule); ApaI GGGCCC at 500, and GGGCCA at
520 and AGGCCC at 540 (one mismatch each). Recorded with
`tools/find_enzyme.sh` and `tools/parse_enzymes.py`.

All 6 screens follow these rules:
- both strands are searched; a site is reported at the top-strand
  position of its leftmost base (its 5' end on the top strand, its
  3' end on the complement strand); complement sites of
  non-palindromic enzymes are marked "C" (true sites) or written in
  parentheses (almost sites); palindromic sites are reported once;
- sites whose cut falls outside a linear molecule are listed;
- sites are sorted by position;
- almost sites are the sites with at most one mismatch, true sites
  included; each shows a dot per matching base and the molecule's
  base where it differs, read on the site's strand in top-strand
  order (a complement site is compared with the reverse complement
  of the recognition sequence, as the header shows: "GGATC
  (GATCC)"). Both lists equal a brute-force search, order included.

More lists (recorded 2026-10-18, same tools; `find_enzyme.sh` with
`DOWNS` set to scroll long lists):
- `enzymes/TENZC.GB`: a circular molecule (LOCUS line "circular"),
  random (seed 4015) with AluI AGCT split across the origin (AG at
  599-600, CT at 1-2) and PleI GAGTC at 300 and GACTC at 400 (plus
  GAGTC at 83 and 372 by chance); PleI's site differs from its
  reverse complement GACTC by one base, so an almost search finds
  sites on both strands at the same position;
- `enzymes/TLONG.GB` (1200 bases) and `enzymes/TLONG3.GB` (3000
  bases): random linear molecules with 63 and 151 almost AluI sites.

All 12 lists (284 sites) follow the rules above, and:
- a circular molecule is searched across the origin: the AluI site
  split across it is reported at 599 (a true site, and an almost
  site "....");
- at the same position, exact sites come before sites with one
  mismatch, then top before complement: PleI almost sites at 300
  (top, exact) then 300 (complement, one mismatch), but 400
  (complement, exact) then 400 (top, one mismatch); at 306, both
  with one mismatch, top first;
- long lists scroll, three sites per line; 63 and 151 sites are all
  listed, without a warning: the "( >50 sites )" and "More than 100
  sites found" messages belong to other screens (probably Find Pair).

The 359-enzyme table cannot be reached from the screens: Find
Enzymes, the F1 = REnz list of the Find Pair form and the F1 = REnz
list of Modify -> Mutagenize (checked 2026-10-18) show the same 89
enzymes of the default menu; INS, DEL and F1-F4 do not change the
menu. No enzyme of that menu has IUPAC letters or N, so how PRIMER.exe
matches them was not observed.

Mutagenize (Modify -> Active Primer: Mutagenize, five downs from the
top of the menu) shows the original and mutated primer with their
translations in three frames, and %GC, Tm, hairpin, dimers and runs of
both. Bases are typed over (overwrite mode; the cursor does not move
by itself). No screen lists restriction sites created by a mutation:
F1 only opens the enzyme list (selecting an enzyme there did not
change the primer; the executable holds the message "Can Only Paste
Enzyme Sequence Into 5'Extension Field"). Not recorded as data: the
key sequence is space, down, enter, enter, the primer, enter, s,
right, down x5, enter, then f1.
