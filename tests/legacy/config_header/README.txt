The first 4 bytes (and bytes 47-54) of pd2.cfg (recorded 2026-10-20)
=====================================================================

Inputs (written with python from ../config/PD2_RUNS42.CFG, an earlier
recording: built-in defaults with runs 4 and 3' G/C runs 2):
  PD2_HEADER_01020304_RUNS42.CFG  bytes 0-3 set to 01 02 03 04
  PD2_TRAILER_05_0C_RUNS42.CFG    bytes 47-54 set to 05 06 ... 0C

Recordings (PRIMER_EXE=... tools/probe_keys.sh WORKDIR NAME "KEYS" FILE,
FILE copied as PD2.CFG; the PD2.CFG left in the run directory is kept
as SAVED_*.CFG):
  SAVED_FROM_HEADER_01020304.CFG, save_from_header_01020304.screens.txt:
    PD2_HEADER_01020304_RUNS42.CFG as PD2.CFG; Modify -> System
    Defaults -> Criteria, runs set to 3 and 3, Quit, y, y (save).
    PAUSE=3, KEYS "space right down down down down down down down enter
    down enter enter enter enter enter enter 3 enter 7 enter 3 3 , esc
    , right right right , enter , , y , , y" (the keys of
    tools/save_run_thresholds.sh). A first run at PAUSE=2.5 lost its
    keys and saved nothing.
  SAVED_FROM_TRAILER_05_0C.CFG, save_from_trailer_05_0c.screens.txt:
    the same with PD2_TRAILER_05_0C_RUNS42.CFG (PAUSE=2.5).
  system_defaults_header_01020304.screens.txt: the three System
    Defaults screens (Defaults, Criteria, Weights) with
    PD2_HEADER_01020304_RUNS42.CFG. PAUSE=2.5, KEYS "space right down
    down down down down down down enter , enter , esc , down enter ,
    esc , down enter , ,".
  SAVED_AFTER_MOLECULE_PRIMER.CFG, save_after_molecule_primer.screens.txt:
    no PD2.CFG; a primer made from ../pairs/TEST600.GB at 51, then runs
    set to 4 and 2 and the defaults saved. PAUSE=3, KEYS "space enter
    enter down down enter , , enter 5 1 enter enter enter enter s ,
    right down down down down down down down enter down enter enter
    enter enter enter enter 3 enter 7 enter 4 2 , esc right right right
    enter , y , y" (file TEST600.GB).

Observations
------------
- Both saved files hold the new run thresholds (offsets 38 and 40: 3
  and 3), and keep the unknown bytes of the file they were loaded from:
  01 02 03 04 at offsets 0-3, 05 to 0C at offsets 47-54. The program
  reads these bytes and writes them back unchanged.
- With 01 02 03 04 in bytes 0-3, the System Defaults screens show the
  values of the other bytes (primer length 20, GC within 5, Tm within
  10, Match -50/+50, Find Pair -100/+100, GC 50-60, Tm 55-80, hairpin
  0.0, dimers 3 and 7, runs 4 and 2, weights 1): no field of these
  screens shows the unknown bytes.
- Saved without a PD2.CFG, after making a primer from a molecule: bytes
  0-3 and 47-54 are zero (SAVED_AFTER_MOLECULE_PRIMER.CFG is byte for
  byte ../config/PD2_RUNS42.CFG).
Not known: whether any action makes PRIMER.exe write non-zero values
there, and what they mean. The menus (recorded while probing: Modify
holds Position, Length, 5'-extension, Criteria, Mutagenize, System
Defaults; Analyze: Active Primer, Primer Pair, Multi Primer;
Molecule-Ops: Find Enzymes, View/Edit, Save/Remove; Quit opens no
menu) show no other settings item; their sub-screens were not all
explored.
