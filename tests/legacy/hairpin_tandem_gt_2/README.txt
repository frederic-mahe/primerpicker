Hairpins with tandem G.T pairs, second series (recorded 2026-10-21)
===================================================================

Observations and fitted values: ../PROBES_20261021_hairpin.md. First
evidence: ../hairpin_tandem_gt/ and ../PROBES_20261021.md, section 6.

117 primers designed to contain hairpins with a tandem G.T step
(5'-TG-3'/3'-GT-5' and 5'-GT-3'/3'-TG-5'), single G.T pairs placed
where the zipping rule matters, and GG/TT steps; typed with Primer/Create
-> Direct Entry, built-in settings (no PD2.CFG), Inspect -> H and D:

  PRIMER_EXE=.../PRIMER.exe tools/inspect_batch.sh LIST RUNDIR PAUSE 6
  (PAUSE 1.5 for the first batch, 2 for the next three, 2.5 for the
  reruns; every primer is kept only from a run that printed "ok"; 20
  runs (18 primers) lost keys and were rerun; one primer never
  completed (three tries) and is left out:
  CGGTGTAGCACACCCCAACAACGCTGTACCC)
  python3 -I tools/parse_inspect.py RUN_DIRS... > hairpin_tandem_gt_2.tsv
  python3 -I tools/hairpin_screens.py RUN_DIRS... > hairpins_screens.txt

(the run directories held screens.txt after tools/strip_splash.py; the
five batches used an ad-hoc script with the same commands as
tools/inspect_batch.sh, which was written from it afterwards and checked
on one primer)

Files:
  hairpin_tandem_gt_2.tsv  parse_inspect.py table (same columns as
                           ../inspect/hairpin_design.tsv)
  designs.tsv              per primer: the design (see below), the shown
                           energy and number of pairs (':' marks), the
                           proposed model (tools/fit_gt_steps.py) and
                           primerpicker at git HEAD ec21a87 (analyze)
  hairpins_screens.txt     the Hairpins screen of each recording (energy
                           and drawing), no splash screen

Design labels (5' arm written outer -> inner; x is the Watson-Crick pair
outside the G.T pairs, y the one inside; outerN/innerN: Watson-Crick
pairs beyond x / y; "sum before" is the range of the step sum allowed by
the recordings available when the primer was designed):
  tandem TG / tandem GT  context x..y: stems ... x T G y ... or ... x G T y
                         ..., at various distances from the loop and from
                         the stem's end, loops 4 to 29 (loop penalties
                         ending in .x0, .x2, .x4, .x5, .x6, .x8 so that the
                         truncation to tenths pins the sums)
  S1   tandem next to the loop (inner end of the stem)
  S2   tandem at the outer end of the stem
  S4   runs of three or four G.T pairs (TGT, GTG, TGTG, GTGT, TGGT, GGTG)
       and two tandems separated by one Watson-Crick pair
  S5   zip test: x, tandem, y, loop of 3 (a stem closed by x with a loop
       of 9 would be more stable than the full stem)
  Z-a  zip test with a single G.T pair: x, G.T, y, loop of 3
  Z-b  x, tandem, y, with a loop of 2 inside y
  Z-c  x, G.T, A.A mismatch, T:A, loop of 3
  Z-d  x, tandem, A.A mismatch, G:C, loop of 3
  GG/TT  x G G y with T T opposite (5'-GG-3'/3'-TT-5')
