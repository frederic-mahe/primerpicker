# Recorded behaviour of the legacy program

Reference observations of *Primer Designer v2.0* (Scientific &
Educational Software, 1990-91), `PRIMER.exe`. The executable itself
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
| 42-46 | five weights, one byte each (order not confirmed) | 1 |
| 47-54 | zero in both files | 0 |

Which of offsets 38 and 40 holds which run threshold is inferred from
the screen order (only offset 38 changed). Navigation notes: in the
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
  rarely appear), and sites are sorted by decreasing Tm (top-strand
  sites first among equal Tm).
Not probed yet: the "exclude Tm <= 0" and "sort by Tm" settings set
to N, and molecules with more sites than one screen holds.


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
