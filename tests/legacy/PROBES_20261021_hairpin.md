# Probes of PRIMER.exe: hairpin stems with tandem G.T pairs (2026-10-21)

Question: how PRIMER.exe (Primer Designer 2.0) scores hairpin stems that hold a
tandem G.T step (G.T then T.G, or T.G then G.T), so that primerpicker's hairpin
model (`core/analysis/legacy_hairpin.hpp`, `core/thermo/freier_1986.hpp`) can be
fixed. First evidence: `PROBES_20261021.md`, section 6, and `hairpin_tandem_gt/`
(three primers, energies -1.0, -6.1 and -5.5 kcal).

Notation. A step is written 5'-ab-3'/3'-cd-5' (a pairs with c, b with d), as in
`freier_1986.hpp`; a step equals its reverse (5'-dc-3'/3'-ba-5'). The two tandem
steps are **TG/GT** (5'-TG-3'/3'-GT-5': T.G then G.T; both arms of the hairpin
read 5'-TG-3') and **GT/TG** (5'-GT-3'/3'-TG-5': G.T then T.G; both arms read
5'-GT-3'). Every step with one G.T pair is one of eight values: L(X), the step
from the Watson-Crick pair X to a G.T pair (G on the 5' strand), or R(Y), the
step from that G.T pair to Y (a T.G pair is a G.T pair read on the other strand):

| table entry (freier_1986.hpp) | role | value at HEAD (hundredths of kcal) |
|---|---|---|
| AG/TT | L(A:T) | 50 |
| CG/GT | L(C:G) | 150 |
| GG/CT | L(G:C) | 130 |
| TA/GT | L(T:A) | 70 |
| GA/TT | R(A:T) | 50 |
| GC/TG | R(C:G) | 190 |
| CT/GG | R(G:C) | 150 |
| AT/TG | R(T:A) | 70 |
| GG/TT | G.T then G.T | 50 |
| GT/TG | tandem G.T then T.G | 50 |
| TG/GT | tandem T.G then G.T | (missing) |

Energies are in hundredths of kcal ("centi"), as the core computes them; PRIMER.exe
shows them truncated toward zero to tenths, so a shown -1.6 means -170 < E <= -160.

## What was run

- 117 primers (16 to 49 bases), typed with Primer/Create -> Direct Entry, built-in
  settings, Inspect -> H (Hairpins) and D, with `tools/inspect_primer.sh` run by
  `tools/inspect_batch.sh` (6 emulators, PAUSE 1.5 to 2.5), on copies of
  `tests/legacy` exported from git HEAD (ec21a87) into `$TMPDIR`. Five batches,
  08:10 to 10:30 (wall time).
- Every primer was designed with a Python copy of primerpicker's model (checked to
  reproduce all 262 earlier recordings), and its prediction written down before the
  recording (`designs.tsv` keeps the design label). Designs were adaptive: after each
  batch, a linear program (`tools/fit_gt_steps.py fit`) gave the range of each step
  sum still allowed by all recordings, and the next primers used loop lengths whose
  penalty ends in .x0, .x2, .x4, .x5, .x6 or .x8 kcal (loops 4-20, 24/29, 23/28,
  13/15, 22/27, 21/26) so that the truncation boundary fell inside that range.
- 20 runs (18 primers) lost keys ("incomplete") and were rerun with a larger PAUSE;
  CGGTGTAGCACACCCCAACAACGCTGTACCC never completed (three tries) and is not recorded.
  No other failure. `git status --short tests/legacy` showed only the new files after
  every batch.

Files (all new):

- `hairpin_tandem_gt_2/hairpin_tandem_gt_2.tsv`: the 117 recordings, columns of
  `inspect/hairpin_design.tsv` (`tools/parse_inspect.py`);
- `hairpin_tandem_gt_2/designs.tsv`: design label, shown energy and pairs, the
  proposed model's energy and pairs, primerpicker HEAD's energy;
- `hairpin_tandem_gt_2/hairpins_screens.txt`: the Hairpins screen of every recording
  (no splash screen);
- `hairpin_tandem_gt_2/README.txt`: commands and design labels;
- `tools/inspect_batch.sh` (parallel `inspect_primer.sh` with splash removal and
  cleanup), `tools/hairpin_screens.py` (keeps the Hairpins screen),
  `tools/fit_gt_steps.py` (`check`: recordings against the model, with
  `--old-zip`/`--no-tandem` for the HEAD rules; `fit`: feasible ranges of the step
  sums, needs scipy).

## Observations

1. **Stems extend through both tandem steps**, in all 16 Watson-Crick contexts
   (x before, y after the tandem) of each orientation, at 0 to 3 pairs from the
   loop and from the stem's outer end, with loops of 3 to 29: 66 designed
   recordings, every drawing shows the full stem (G.T pairs marked `:`).
2. **A tandem next to the loop is not part of the stem** (S1, 6 recordings: x=A, C,
   T, both orientations): the stem ends at x and the two G.T pairs join the loop
   (for x=A: -3.9 shown; with the tandem it would be -5.0). Same rule as for a single
   G.T pair (README: "never at either end").
3. **A tandem at the outer end of the stem is dropped** (S2, 7 recordings): the stem
   starts at y (for y=T: -3.8 shown, 4 pairs; with the tandem: -5.1).
4. **Stems are zipped inwards through G.T pairs.** With a stem x, tandem, y and a
   loop of 3 inside y (S5, 10 recordings), HEAD's rule (zip over Watson-Crick pairs
   only) also scores the hairpin closed by x with a loop of 9, which is more stable,
   and shows it (4 pairs, e.g. -3.9). PRIMER.exe shows the full 7-pair stem with the
   loop of 3 (-2.5 for AGCCATGAAAATTGTGGCA). The same holds for a **single** G.T pair
   (Z-a, 8 recordings: x, G.T, y, loop 3: PRIMER.exe shows 6 pairs, -1.7 for
   AGCCAGAAAATTTGGCA, where HEAD shows 4 pairs, -4.0). None of the 262 earlier
   recordings distinguishes the two rules.
5. Limits of the zipping: an inner Watson-Crick pair that would leave a loop of 2 does
   not count (Z-b, 4 recordings: the stem closes outside the tandem, loop 8, -4.0 /
   -3.9), and zipping stops at a mismatch (Z-c, Z-d: x, G.T or tandem, A.A, then a
   Watson-Crick pair: the stem closed by x is shown, -3.9 / -3.5). These four agree
   with HEAD's rule too; they rule out a "zip to any inner pair" variant.
6. **Runs of G.T pairs are additive nearest-neighbour steps** (S4, 9 recordings:
   TGT, GTG, TGTG, GTGT, TGGT, GGTG runs, two tandems separated by one Watson-Crick
   pair): all equal the sums of the steps (TG/GT, GT/TG, GG/TT) fitted elsewhere.
7. GG/TT (5 recordings with phases .x0 to .x6) is pinned at 0.50 (HEAD's value).

## Fitted values

Each recording gives E = loop penalty - (Watson-Crick stacks) - (G.T steps), with
-10 < E - 10 x shown <= 0 (hundredths). Since the eight single-G.T steps are only
known through sums L(X) + R(Y) (adding c to every L and subtracting it from every R
changes no prediction, and then GT/TG changes by +2c and TG/GT by -2c), the
recordings determine sums, not single values. `tools/fit_gt_steps.py fit` over the
262 earlier recordings, the three of `hairpin_tandem_gt` and `find_pair_trun` (1867),
and the 117 new ones (14653 inequalities, all satisfiable) gives these ranges (in
hundredths; the upper ends are bounds of a continuous relaxation):

| sum | range | with HEAD's table + TG/GT = 60 | range before these probes |
|---|---|---|---|
| R(A:T) + TG/GT + R(A:T) | [160, 161] | 160 | [151, 169] |
| R(A:T) + TG/GT + R(C:G) | [300, 301] | 300 | [295, 309] |
| R(A:T) + TG/GT + R(G:C) | [260, 261] | 260 | [255.5, 269] |
| R(A:T) + TG/GT + R(T:A) | [180, 181] | 180 | [171.5, 198] |
| R(C:G) + TG/GT + R(C:G) | [440, 442] | 440 | [436, 458] |
| R(C:G) + TG/GT + R(G:C) | [400, 402] | 400 | [400, 409] |
| R(C:G) + TG/GT + R(T:A) | [320, 321] | 320 | [316, 338] |
| R(G:C) + TG/GT + R(G:C) | [360, 362] | 360 | [360, 369] |
| R(G:C) + TG/GT + R(T:A) | [280, 281] | 280 | [276, 298] |
| R(T:A) + TG/GT + R(T:A) | [200, 202] | 200 | [192, 227] |
| L(A:T) + GT/TG + L(A:T) | [150, 151] | 150 | [144, 159] |
| L(A:T) + GT/TG + L(C:G) | [250, 251] | 250 | [247, 259] |
| L(A:T) + GT/TG + L(G:C) | [230, 231] | 230 | [230, 239] |
| L(A:T) + GT/TG + L(T:A) | [170, 171] | 170 | [166, 181.5] |
| L(C:G) + GT/TG + L(C:G) | [350, 351] | 350 | [350, 359] |
| L(C:G) + GT/TG + L(G:C) | [330, 331] | 330 | [330, 342] |
| L(C:G) + GT/TG + L(T:A) | [270, 271] | 270 | [270, 279] |
| L(G:C) + GT/TG + L(G:C) | [310, 312] | 310 | [310, 328] |
| L(G:C) + GT/TG + L(T:A) | [250, 251] | 250 | [250, 263.5] |
| L(T:A) + GT/TG + L(T:A) | [190, 191] | 190 | [186, 204] |
| GT/TG + TG/GT (no gauge freedom) | [110, 112] | 110 | [102, 128] |
| GG/TT | [49.5, 51] | 50 | [48, 59] |

The single-G.T sums L(X) + R(Y) are pinned as well (all 16 within [v, v + 2], v =
HEAD's table: L(A:T)+R(A:T) = 100 exactly, ..., L(T:A)+R(T:A) = [140, 141.5]),
because the tandem designs also cross them.

So, **keeping HEAD's ten values, TG/GT = 60 hundredths (0.60 kcal)**: every TG/GT
sum gives TG/GT in [60, 61] (for instance R(A:T) + TG/GT + R(A:T) in [160, 161]
with R(A:T) = 50), and GT/TG = 50 (HEAD's value) is confirmed in all ten contexts.
Every sum is at the lower end of its range at a multiple of 10, which is what a
table with one decimal gives. The 10-context consistency of each tandem step shows
that PRIMER.exe treats the tandem as two nearest-neighbour steps plus one tandem
value, with no context dependence beyond nearest neighbours.

Which recordings pin each sum: Appendix (each sum has 2 to 6 recordings of this
series with different loop phases; the first evidence adds 1867 and
ACCTGGCAAAAGCTGGGTTA to R(C:G) + TG/GT + R(G:C) and ACCGTGGAAAACCGTGGTTA to
L(C:G) + GT/TG + L(C:G)).

## Zipping rule (proposed)

A hairpin is scored for each Watson-Crick pair (i, j) closing a loop of at least 3,
unless, going inwards from it over Watson-Crick and G.T pairs (stopping at the first
other pair), one reaches a Watson-Crick pair that still closes a loop of at least 3:
then the stem is scored there instead. The stem is extended outwards over
Watson-Crick and G.T pairs (all eleven steps now exist), and G.T pairs at its outer
end are dropped. G.T pairs next to the loop stay in the loop. (HEAD: inwards over
Watson-Crick pairs only.)

## Every recording against the proposed model

`tools/fit_gt_steps.py check` (TG/GT = 60, zipping through G.T pairs):

| recordings | proposed | HEAD rules (`--old-zip --no-tandem`) | tandem only (`--old-zip`) | zip only (`--no-tandem`) |
|---|---|---|---|---|
| 117 of `hairpin_tandem_gt_2` | 0 misfits | 58 | 18 | 45 |
| 262 of `inspect/` (inspect, hairpin_design, inspect_2) | 0 | 0 | 0 | 0 |
| 3 first evidence (1867 and `hairpin_tandem_gt/`) | 0 | 2 | 0 | 2 |

A misfit is a different shown energy or a different number of pairs. The same model in
C++ (the change below, applied to a copy of the repository in `$TMPDIR`, not to the
tree) passes all 256 ctest tests of the debug preset, `legacy hairpins match the
recorded primers` with the 117 new recordings added (380 assertions), and `hairpin
drawings match every recorded drawing` with them added (285 drawings: the 5' line,
the marks and the 3' line of every new recording are reproduced). The patched CLI
gives -1.0, -6.1 and -5.5 for the three first-evidence primers.

**Misfits: none.**

## Not settled

- Single step values: only the sums above are determined (gauge freedom of the eight
  single-G.T steps). TG/GT = 0.60 holds for HEAD's choice of the ten other values.
  Whether these numbers are the printed Freier et al. 1986 G.U values was not
  checked: the paper is not in the repository.
- Hundredths: each sum is pinned to [v, v + 1] or [v, v + 2] hundredths, not to
  exactly v (no loop penalty ends in .x1 or .x9 below 40 bases of loop).
- Runs of three or four G.T pairs and two tandems: one recording each (consistent
  with additivity, not each pinned twice).
- Zipping through a G.T pair when the next Watson-Crick pair is two or more G.T pairs
  further in was recorded with tandems (S5) and single pairs (Z-a); longer G.T runs
  before the inner Watson-Crick pair were not tested.
- Hairpins of the Multi Primer screen and Find Pair's hairpin criterion were not
  re-recorded; they are assumed to use the same function (PROBES_20261021.md,
  section 6, showed Find Pair finds the 1867 hairpin).

## Draft for README.md (replaces the last two items of the hairpin model)

> - stems are zipped inwards through Watson-Crick and G.T pairs, to the innermost
>   Watson-Crick pair that keeps a loop of at least 3 bases (zipping stops at a
>   mismatch); only these maximal stems are scored (`ATAT...` repeats get a loop of
>   4, not 6; `hairpin_tandem_gt_2/`, designs S5 and Z-a: a stem closed outside a
>   G.T pair with a larger loop is not scored when a Watson-Crick pair further in
>   closes a loop of 3);
> - ...
> - stems extend through tandem G.T steps (recorded 2026-10-21,
>   `hairpin_tandem_gt_2/`, `PROBES_20261021_hairpin.md`): 5'-TG-3'/3'-GT-5' has
>   its own value (0.60 kcal with the other G.T values of `freier_1986.hpp`), and
>   5'-GT-3'/3'-TG-5' (0.50) and GG/TT (0.50) are confirmed; each tandem step was
>   recorded in the 16 Watson-Crick contexts, all consistent with nearest-neighbour
>   additivity; as for single G.T pairs, tandem G.T pairs at the outer end of a stem
>   or next to the loop are not part of the stem.

## Proposed change for freier_1986.hpp (values only)

Add the tandem step to `freier_1986_wobble_table`:

    WobbleStack{Step{'T', 'G', 'G', 'T'}, 60},

Source: fitted on PRIMER.exe recordings (`tests/legacy/hairpin_tandem_gt_2/`,
`PROBES_20261021_hairpin.md`); with the ten other values unchanged, all TG/GT sums
give 60 or 61 hundredths, the shown energies all fit with 60. The comment above the
table ("The tandem step T.G over G.T is not in the table: stems do not extend
through it") and the zipping comments in `legacy_hairpin.hpp` ("a stem that could zip
one more pair inwards is scored there", "(over Watson-Crick pairs only)") then no
longer describe the behaviour: changing them is left to a human. The zipping change
for `legacy_hairpin.hpp` that was checked in the copy:

```diff
+// true if a Watson-Crick pair further in, reached through Watson-Crick
+// and G.T pairs, still closes a loop of at least three bases
+constexpr auto zips_inward(std::span<sequence::Nucleotide const> const primer, std::size_t const inner_5prime,
+                           std::size_t const inner_3prime) noexcept -> bool {
+    for (auto depth = std::size_t{1}; inner_5prime + depth < inner_3prime - depth; ++depth) {
+        auto const top = primer[inner_5prime + depth];
+        auto const bottom = primer[inner_3prime - depth];
+        if (inner_3prime - depth - (inner_5prime + depth) - 1 < thermo::freier_1986_min_hairpin_loop) {
+            return false;
+        }
+        if (pairs(top, bottom)) {
+            return true;
+        }
+        if (not is_wobble(top, bottom)) {
+            return false;
+        }
+    }
+    return false;
+}
 ...
             // a stem that could zip one more pair inwards is scored there
-            auto const loop_length = inner_3prime - inner_5prime - 1;
-            if (loop_length >= min_loop + 2 and detail::pairs(primer[inner_5prime + 1], primer[inner_3prime - 1])) {
+            if (detail::zips_inward(primer, inner_5prime, inner_3prime)) {
                 continue;
             }
```

(The two unit tests that read the recordings would need
`hairpin_tandem_gt_2/hairpin_tandem_gt_2.tsv` added to their file lists, and the
drawing test's count 168 -> 285; existing tests were not changed in the tree.)

## What went wrong

- 20 runs lost keys (mostly primers of 35 to 49 bases, or DOSBox-X stopping early
  while typing); all but one primer were recorded on a rerun.
- My first linear-program fit after the zip recordings was infeasible because it
  still enumerated stems with HEAD's zipping rule; with the proposed rule it is
  feasible (this is how observation 4 was confirmed as a rule, not a value change).
- `tools/inspect_batch.sh` was written after the batches from the ad-hoc script
  used in `$TMPDIR` (same commands); it was checked on one primer only.

## Appendix: recordings per step sum

Shown energy and loop length; ":" = the shown stem has no G.T step (S1, S2, Z-b, Z-c,
Z-d, where the G.T pairs are outside the stem).

```
: 19
    AGCCGTAAAAGTGGCA (-2.2, loop 8)
    AGCCTGAAAATGGGCA (-2.2, loop 8)
    AGTGCCGAAACGGCGTA (-.9, loop 3)
    ATGGCCGAAACGGCTGA (-.9, loop 3)
    AGCCAGTCAAGGTTGGCA (-4.0, loop 8)
    AGCCATGCAAGTGTGGCA (-4.0, loop 8)
    AGGCTGTGCACGTAGCCA (-3.9, loop 8)
    AGGCTTGGCACTGAGCCA (-3.9, loop 8)
    AGTCGCGAAAACGCGGTA (-1.5, loop 4)
    AGCCAGATCCCAATTGGCA (-3.9, loop 9)
    AGCCAGTAAAAAGTTGGCA (-3.9, loop 9)
    AGCCATGAAAAATGTGGCA (-3.9, loop 9)
    AGTACGGAAAAACCGTGTA (-2.6, loop 5)
    ATGACGGAAAAACCGTTGA (-2.6, loop 5)
    AGGCTGTAAAAAAGTAGCCA (-3.7, loop 10)
    AGGCTTGAAAAAATGAGCCA (-3.7, loop 10)
    AGTTGGCAAAAAAGCCAGTA (-3.8, loop 6)
    ATGTGGCAAAAAAGCCATGA (-3.8, loop 6)
    AGCCATGAGCCCAATGTGGCA (-3.5, loop 11)
GT/TG + L(T:A) + L(T:A): 3
    CGTGTAGACCCAACCAACCCTGTACA (-.8, loop 12)
    CGCTGTACACCACCCAAACAACCCACAAACACGTGTAGCC (-1.2, loop 24)
    CCCTGTAGCCACACAAACACCACAACACCCACCGCTGTAGGA (-3.7, loop 24)
L(A:T) + GT/TG + L(T:A): 6
    AGCCAGTAAAATGTTGGCA (-2.4, loop 3)
    AGGCTGTTCCCAGTAGCCA (-2.3, loop 3)
    ATGTTCGCCCCCCGAGTAA (-1.6, loop 5)
    CCGGAGTAGAAACACACAACACCCCCTCTGTTCCGC (-6.8, loop 16)
    CGTGTTCGCACCCCCCCCAAACCACACCACCCGAGTACC (-.4, loop 23)
    CGAGTAAGGAAAAAACACCCCCACCCCCCCCACACCCCCTTGTTCC (-.9, loop 28)
L(A:T) + L(A:T) + GT/TG: 4
    ACGGAGTTAAAAGTTCCGA (-1.3, loop 3)
    AGAGTTCAACCACCACACCAGAGTTCA (-.8, loop 13)
    CCCAGTTCCACACAAACAAAAACAAAAAACCAACCGGAGTTGGC (-3.2, loop 26)
    CGCAGTTCCCCCAAACACCCCCAAACACAACCACACAAGGAGTTGCA (-3.1, loop 29)
L(A:T) + L(C:G) + GT/TG: 4
    CGCGTTCAACCAGAGTGCA (-3.8, loop 5)
    CGGTAGTGCCACCCCAAAAAAACACAACACGGCGTTACCA (-7.8, loop 20)
    ATGCAGTGCCACACCAAACACACACCCACCACGCGTTGCAC (-5.2, loop 23)
    AGGCGTTGCAACCAAACCAAACCCACCAAAACCCCCCGCAGTGCCA (-5.4, loop 28)
L(A:T) + L(G:C) + GT/TG: 4
    ACCGGGTTGCAAAACACACACACGCAGTCCGGA (-10.0, loop 13)
    CCGAGTCCACCCCACCACAAAAAGTGGGTTCGA (-8.1, loop 13)
    AGGGGTTGGCCACACAAAACACACAAAAACCACCGCCAGTCCCA (-8.3, loop 24)
    ACGAGTCGGCACACAAACCAACACACACCCACACAAACCCGGTTCGC (-2.7, loop 29)
L(A:T) + R(A:T): 2
    AGCCAGAAAATTTGGCA (-1.7, loop 3)
    AGGCTTTCCCAGAGCCA (-1.6, loop 3)
L(A:T) + R(G:C) + GG/TT: 1
    ACGAGGGCAACCCCCACAAACACCAAACCACCGCTTTCGA (-2.3, loop 24)
L(A:T) + R(T:A): 2
    ACGGAGTAAAATTCCGA (-1.0, loop 3)
    ACGGATTAAAAGTCCGA (-1.0, loop 3)
L(C:G) + GT/TG + L(T:A): 5
    AGCGCGTACAATGTGCGCA (-4.1, loop 3)
    CGCGTAGACACACAACCCCACTGTGCA (-2.5, loop 13)
    ACCTGTGGCAAAACCCCCAAACACGCCGTAGGC (-7.7, loop 15)
    ATGTGGGCAACAACAACCCCAAAAAACAACCAGCCCGTAC (-4.0, loop 24)
    AGGGCGTACCACCCCAACAAAAACAAACAAAAACCTGGTGTGCCCA (-10.8, loop 24)
L(C:G) + L(C:G) + GT/TG: 2
    ACGTGCGACCAAACAAAAAACGCGTGA (-3.6, loop 13)
    CGCCGTGTCCAAACCACCCAAACCAACCCACAAACAACGACGTGGCA (-5.4, loop 29)
L(C:G) + L(C:G) + GT/TG + GT/TG + TG/GT: 1
    AGGCGTGTGCAAAAGCGTGTGCCA (-8.4, loop 4)
L(C:G) + L(C:G) + L(G:C) + L(G:C) + GT/TG + GT/TG: 1
    ACCGTCGTCCAAAAGGGTGGTGGA (-6.5, loop 4)
L(C:G) + L(C:G) + R(G:C) + R(C:G) + GT/TG + TG/GT: 1
    AGCTGCGTGCAAAAGCGTGTGGCA (-8.4, loop 4)
L(C:G) + L(G:C) + GT/TG: 4
    CAGGTGATCCAAAACACAAACAACCAAAAAATCGTCTA (-.7, loop 22)
    CGGTGGCGACAAACAACCCCAACCCCAACAAACGCCGTCA (-3.7, loop 24)
    ACGCGTCGGCCAAACCAAACAACACCCCAAACCCCGGTGCGA (-5.7, loop 24)
    CCGCGTCGGGCCCACACCAACAAAAAAAACACCCCGGTGCGC (-9.0, loop 22)
L(C:G) + R(A:T): 1
    AGCGCGACAATTGCGCA (-3.4, loop 3)
L(C:G) + R(A:T) + GG/TT: 1
    ACGCGGAGCAAAACCAAACCAAAAAAACCCCCTTTGCGC (-1.9, loop 23)
L(C:G) + R(C:G) + GG/TT: 1
    CCCGGCGCACACCACCCCACACCCGTTGGA (-2.7, loop 16)
L(C:G) + R(C:G) + GG/TT + GT/TG + TG/GT: 1
    ACGCGGTGCCAAAAGGTGTTGCGA (-7.4, loop 4)
L(C:G) + R(C:G) + GT/TG + TG/GT: 1
    AGCGTGCCAAAAGGTGTGCA (-4.9, loop 4)
L(C:G) + R(G:C) + GG/TT + GT/TG + TG/GT: 1
    AGCCTGGTGCAAAAGCGTTGGGCA (-8.4, loop 4)
L(G:C) + GT/TG + L(T:A): 5
    ACCGGGTACCATGTCCGGA (-2.9, loop 3)
    CCTGTCGCCCCCACCAAAACAACCCCACACGCGGTAGA (-2.1, loop 22)
    ACCCTGTCGCCCACCCCACCCAAAAAACCAAAACGGTAGGGC (-4.1, loop 24)
    CCCGGTACGAAACACACAACCAAACACAAAAACAAACGTGTCGGA (-3.1, loop 27)
    AGAGGTATGCCACAAACAACAACAAAACAAACCCCAACCATGTCTCA (-.4, loop 29)
L(G:C) + L(G:C) + GT/TG: 2
    CGGTCGGACAACCCGGTCA (-3.6, loop 5)
    CGCGGTCGACAAAACACAACACACGGTCGCC (-4.6, loop 15)
R(A:T) + L(G:C) + GG/TT: 1
    ACCGGGACAACACCACCAACCCACCCACCCGTTTCGGA (-1.8, loop 22)
R(A:T) + R(A:T) + TG/GT: 2
    ATTGACTCACCCAACACAAGAGTTGAA (-3.1, loop 11)
    ACTTGAAGGCACCACCACCCCCACCCACCTTTGAGC (-2.1, loop 18)
R(A:T) + R(C:G) + TG/GT: 5
    ACCGGTGACCATTGCCGGA (-3.4, loop 3)
    CATTGCGGCACCCACCGTGATC (-4.5, loop 6)
    CCGGGTGAGAAACACACAACACCCCCTCTTGCCCGC (-8.7, loop 16)
    AGGCTTGCCCCACACAACAACACCACAAACCAACGGGTGAGCCA (-8.9, loop 24)
    AGGGTGAGCACAAACAACCCCACACAACCACACACCCCGCTTGCCCC (-5.1, loop 29)
R(C:G) + R(C:G) + TG/GT: 2
    ACGTGCCAACCCAGTGCGA (-2.3, loop 7)
    AGCCGTGCCAAAACCACCCACAACCAAACCCCCCGTGCGGCA (-4.5, loop 26)
R(G:C) + L(G:C) + GG/TT: 1
    ATGGGGACCAACACACCACCTCTTCAA (-2.1, loop 13)
R(G:C) + L(G:C) + GT/TG + TG/GT: 1
    AGCTGTCGAAAACGGTGGCA (-3.4, loop 4)
R(G:C) + L(T:A): 1
    AGCGCTACAATGGCGCA (-3.6, loop 3)
R(G:C) + R(A:T) + TG/GT: 5
    AGCGCTGACAATTGGCGCA (-4.0, loop 3)
    CGTTGGCCCACCAAACTGACC (-.5, loop 9)
    ATGCTGAGTCACCCAAACAACCCACTTGGCAC (-6.0, loop 14)
    AGGCCTGACAGAACCAAAAACACACAAAAAAACACCTGTTGGGCCA (-9.5, loop 24)
    CGATTGGGCGACACCACCCAACCCCACCCAACAACCCCCCGCCTGATCC (-5.3, loop 29)
R(G:C) + R(C:G) + GT/TG + TG/GT + TG/GT: 1
    AGCCTGTGCCAAAAGGTGTGGGCA (-8.4, loop 4)
R(G:C) + R(C:G) + TG/GT: 4
    CCCTGCCAAACCAACCCCCGTGGGC (-1.6, loop 13)
    AGGTGGGGCACACCAAACCCCTGCCA (-8.4, loop 10)
    CCCTGCCTCAACCAACCAAGGTGGGC (-7.2, loop 10)
    ATCGGTGGCAAAAACAACACCCTGCCGAC (-5.9, loop 13)
R(G:C) + R(G:C) + R(G:C) + R(C:G) + TG/GT + TG/GT: 1
    AGCTGCTGGCAAAAGCTGGTGGCA (-8.5, loop 4)
R(G:C) + R(G:C) + TG/GT: 2
    CCACTGGGCCCCCCCAAAAGGCCTGGTGC (-12.5, loop 9)
    CGCTGGGCCAACACCAAACCCCCGCCTGGCA (-7.4, loop 15)
R(T:A) + L(T:A): 2
    AGCCATAAAATGTGGCA (-2.1, loop 3)
    AGGCTGTCCCATAGCCA (-2.0, loop 3)
R(T:A) + R(A:T) + TG/GT: 6
    AGCCATGAAAATTGTGGCA (-2.5, loop 3)
    AGGCTTGTCCCATGAGCCA (-2.4, loop 3)
    ACTTGTCCCCACACCCACCCGATGAGA (-.5, loop 13)
    CGCATGAGGAAAACCCCAAACAACCCTTGTGCC (-5.7, loop 15)
    ACCTTGTGCAACACCAACACACAAACCACAACCACACAGCATGAGGA (-2.8, loop 29)
    CGGATGACTCACCAAAAAAAAACAAAACCCCACAAACAAGTTGTCCA (-2.0, loop 29)
R(T:A) + R(C:G) + TG/GT: 4
    AATGCGCAAACCCAAAAAACCAACCACCACCGCGTGTC (-.7, loop 24)
    AGGATGCCCCCCCCACCCCACCAACCAAACACACGGTGTCCA (-3.1, loop 26)
    ACCGGTGTAGAACCCCCCCAAAACCACACAAAAACTATGCCGGA (-5.9, loop 24)
    CGGTGTCCGCACACCACACAAAAAACACCCCCACCCGGATGCCC (-5.1, loop 26)
R(T:A) + R(G:C) + TG/GT: 4
    CGCATGGCAACCCACCCAACCCCACCACAGCTGTGCC (-4.1, loop 21)
    ACCTGTCGCACAACCCCAAAACCAAAAAACACCCCGATGGGC (-1.8, loop 26)
    AGGCTGTTACCCAAACAACCCACACCACCACCCCAACCTAATGGCCA (-2.3, loop 29)
    CGCATGGCGAACACACAAACCCACCCACCCCACAAAAACGCTGTGCA (-4.6, loop 29)
R(T:A) + R(T:A) + TG/GT: 3
    ACGGATGTAAAATGTCCGA (-1.8, loop 3)
    ACCGATGTGGAAAAACACCATGTCGGC (-9.8, loop 7)
    CCCGATGTACCACACACCCACCACAAACACCATGTCGGC (-1.5, loop 23)
```
