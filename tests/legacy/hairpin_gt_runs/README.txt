Hairpins with runs of three or four G.T pairs (recorded 2026-10-22)
=====================================================================

32 primers (20 to 23 bases) typed with Primer/Create -> Direct Entry,
built-in settings, then Inspect, H and D, recorded with
  PRIMER_EXE=... tools/inspect_batch.sh designs.tsv RUNDIR 2 4
(PAUSE 2, 4 emulators; the 4 runs that lost keys again with PAUSE 2.5
and 2 emulators). Each primer is a hairpin T + arm + loop (AAAA or
AAAAA) + arm' + A whose stem holds a run of 3 or 4 G.T pairs (GGG, TTT,
GTG, TGT, GGT, TTG, GTT, TGG, GGGG, TTTT, GTGT, TGTG) between C:G
pairs, in the middle of the stem, next to the loop, or at its outer
end; arm' pairs with the arm, T opposite each G and G opposite each T
of the run.

Files:
- designs.tsv: the primer, its design, and primerpicker's hairpin
  energy computed before the recording (inputs and predictions);
- hairpin_gt_runs.tsv: tools/parse_inspect.py of the 32 runs (columns
  of ../inspect/hairpin_design.tsv);
- hairpins_screens.txt: the Hairpins screen of each run
  (tools/hairpin_screens.py on the runs' screens, splash screens
  removed by tools/strip_splash.py).
All 32 shown energies equal the predictions.
