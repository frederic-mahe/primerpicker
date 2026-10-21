Find Pair counts on search/TRUN.GB, and the 1000-pair cap (recorded 2026-10-21)
===============================================================================

Observations: ../PROBES_20261021.md, sections 6 and 2.

Inputs (not recordings), written for these probes:
  PD2_RUNS4_ADJ5.CFG, PD2_RUNS4_ADJ6.CFG, PD2_RUNS4_ADJ8.CFG:
    python3 -I tools/make_pd2_cfg.py --criteria 50,60,55,80,0,3,A,4,3
    (A = 5, 6, 8; with A = 7 the tool writes ../config/PD2_RUNS4.CFG)
  PD2_CRITERIA_W10131.CFG: ../config/PD2_CRITERIA.CFG with bytes 42-46 set
    to 1 0 1 3 1 (weights GC 1, Tm 0, 3' dimers 3, any dimers 1, false
    priming 1)

Recorded with
  CFG=PD2_RUNS4_ADJ6.CFG tools/print_report.sh WORKDIR find_pair search/TRUN.GB 1280 1380
  CFG=PD2_RUNS4_ADJ8.CFG tools/print_report.sh WORKDIR find_pair search/TRUN.GB 1280 1380
  CFG=PD2_RUNS4_ADJ5.CFG tools/print_report.sh WORKDIR find_pair search/TRUN.GB 880 980
  CFG=config/PD2_RUNS4.CFG PAUSE=4 tools/find_pair_limits.sh WORKDIR search/TRUN.GB 1680 100 B OFFSET 600
    (B OFFSET = 1784 0, 1821 0, 1867 0, 1867 3)
  CFG=config/PD2_CRITERIA.CFG PAUSE=4 tools/find_pair_limits.sh WORKDIR pairs/TEST600.GB 71 0 422 0 600
    (and 73 0 422 0: a control)
  CFG=PD2_CRITERIA_W10131.CFG tools/print_report.sh WORKDIR find_pair pairs/TEST600.GB 260 400
  CFG=config/PD2_RUNS4.CFG tools/print_report.sh WORKDIR inspect ATCTGTGGGACCGCCTGCTT
  PAUSE=3 tools/inspect_primer.sh WORKDIR ATCTGTGGGACCGCCTGCTT   (built-in settings)
(PAUSE 3 unless noted). Files: NAME.screens.txt (splash screens removed by
tools/strip_splash.py) and NAME.PRN (the report), NAME giving the molecule,
targets (or targets and offsets) and settings:
  TRUN_1280_1380_RUNS4_ADJ6     A 22, B 20, 173 pairs; .PRN
  TRUN_1280_1380_RUNS4_ADJ8     A 22, B 27, 300 pairs (the report was not
                                written: keys lost after "Print best")
  TRUN_880_980_RUNS4_ADJ5       A 19, B 29, 160 pairs; .PRN
  TRUN_1680_100_1784_0_RUNS4    A 13, B 1, 1 pair
  TRUN_1680_100_1821_0_RUNS4    A 13, B 1, 1 pair
  TRUN_1680_100_1867_0_RUNS4    "No Matching Primers Found"
  TRUN_1680_100_1867_3_RUNS4    A 13, B 1, 6 pairs (primerpicker: 2 and 7)
  TEST600_71_0_422_0_CRITERIA   A 1, B 1, 1 pair
  TEST600_73_0_422_0_CRITERIA   A 1, B 1, 1 pair
  TEST600_260_400_CRITERIA_W10131  A 37, B 45, "Pairs: 1000"; .PRN (best 20;
                                the 1001st pair generated, 260/450, is row 7)
  INSPECT_ATCTGTGGGACCGCCTGCTT_RUNS4  One Primer report of primer B 1867:
                                "No Hairpins ... NO"
  INSPECT_SCREENS_ATCTGTGGGACCGCCTGCTT  Primer Analysis, Hairpins (energy -1.0,
                                a stem with a tandem T.G/G.T step) and Dimers
                                screens of the same primer
