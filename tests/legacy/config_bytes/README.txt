pd2.cfg bytes 0-3 and 47-54 (recorded 2026-10-21)
==================================================

Observations: ../PROBES_20261021.md, section 5. Earlier probe: ../config_header/.

Recording (a save):
  SAVED_DEFAULTS_NUMERIC2.CFG, save_defaults_numeric2.screens.txt: no PD2.CFG at
  the start; Modify -> System Defaults -> Defaults changed (length 22, GC 71,
  Tm 2, Match Active Primer 33 / 44, Find Pair 66 / 77), quit, save:
    PRIMER_EXE=... PAUSE=3.5 tools/probe_keys.sh WORKDIR save_defaults_numeric2 \
      "space right down down down down down down down enter enter 2 2 enter enter
       7 1 2 enter 3 3 enter 4 4 enter 6 6 enter 7 7 enter , esc , right right
       right , enter , , y , , y"
  (an earlier run with "7 enter" for GC, and one meant to set exclude to N,
  saved nothing: keys went to other fields; not kept).

Inputs (not recordings), written with python from tools/make_pd2_cfg.py
(built-in defaults, PD2_BASE.CFG) or from ../config/PD2_RELAXED.CFG
(PD2_RELAXED_*.CFG), with bytes 0-3 / 47-54 set to:
  *_F1_D1    float32 1.0 (00 00 80 3F) / float64 1.0
  *_F50_D50  float32 50.0 (00 00 48 42) / float64 50.0
  *_FF       FF FF FF FF / FF x 8
  *_I10      int32 10 (0A 00 00 00) / 0A x 8

Recordings with each input as PD2.CFG (PAUSE 3), V being the input's name
without PD2_:
  V_inspect_CTCAAAGGCTCAGACATCAG.PRN   CFG=... tools/print_report.sh WORKDIR inspect CTCAAAGGCTCAGACATCAG
  V_find_pair_TEST600_150_400.PRN      CFG=... tools/print_report.sh WORKDIR find_pair pairs/TEST600.GB 150 400
  V_homology_TREP_51_u.PRN             CFG=... tools/print_report.sh WORKDIR homology homology/TREP.GB 51 u
  V_pair_product_TEST600_51_400c.screens.txt
      tools/probe_keys.sh WORKDIR pair_product "space enter enter down down down
      enter , , enter 5 1 enter enter enter enter s , enter enter 4 0 0 enter y
      enter enter enter s , right right down enter 1 enter 2 enter a , , ,"
      pairs/TEST600.GB PD2.CFG   (the keys of tools/pair_product.sh, one more
      down for PD2.CFG; the Analyze Amplified Product screen)
Files missing from this list were not recorded (keys lost).
