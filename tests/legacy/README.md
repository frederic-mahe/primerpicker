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
