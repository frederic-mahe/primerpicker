Modify -> Length (recorded 2026-10-21)
======================================

Observations: ../PROBES_20261021_picking.md, "Modify -> Length".

A primer made with From Molecule on ../pairs/TEST600.GB, saved with S
(entry 1 of the Primer List, active), then Modify -> Length (Modify,
two downs, Enter): the CHANGE PRIMER LENGTH form ("Primer length (incl
5' ext): 20", type over the number, Enter). Recorded with
  PRIMER_EXE=... PAUSE=3 tools/probe_picking.sh WORKDIR NAME "$BASE KEYS" ../pairs/TEST600.GB
  BASE="space enter enter down down enter , , enter"
(PAUSE 3.5 for length_575_30 and length_51_1, whose first runs froze):

  NAME            KEYS
  length_51       "5 1 enter enter enter enter , s , , right down down enter , , , , ,"
  length_51_25    "5 1 enter enter enter enter , s , , right down down enter , , 2 5 enter , , , , ,"
  length_51_0     "5 1 enter enter enter enter , s , , right down down enter , , 0 enter , , , , ,"
  length_51_1     "5 1 enter enter enter enter , s , , right down down enter , , 1 enter , , , , ,"
  length_575_30   "5 7 5 enter enter enter enter , s , , right down down enter , , 3 0 enter , , , , ,"
  length_400c_25  "4 0 0 enter y enter enter , , s , , right down down enter , , 2 5 enter , , , , ,"

Files: NAME.screens.txt (tools/screen_text.py, splash screens removed by
tools/strip_splash.py). Not kept: the first runs of length_575_30 and
length_51_1 (froze before the form).
