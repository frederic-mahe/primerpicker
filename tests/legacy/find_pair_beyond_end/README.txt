Find Pair targets near and beyond the molecule's end (recorded 2026-10-21)
==========================================================================

Molecule ../pairs/TEST600.GB (600 bases). Settings: ../config/PD2_RELAXED.CFG
copied as PD2.CFG (files named *_RELAXED), or none (built-in settings, files
named *_builtin). Observations and rule: ../PROBES_20261021.md, section 1.

Recorded with
  PRIMER_EXE=... PAUSE=3 tools/probe_keys.sh WORKDIR NAME "KEYS" TEST600.GB [PD2.CFG]
KEYS, for target A TA, offset OA, target B TB, offset OB (digits one by one;
a 2-digit offset needs enter, a 3-digit one moves on by itself), for example
TEST600_150_40_590_40_RELAXED:
  space down down down enter enter down down down enter , , enter 1 5 0 enter
  4 0 enter enter enter 5 9 0 enter 4 0 enter enter enter enter enter , , , ,
  space , , p , enter , down down enter , r e p o r t enter , ,
(without PD2.CFG: "down down" before the molecule's enter).
TEST600_600_40_600_40_RELAXED was recorded with the same keys up to
"space , ," (no printing keys).

Files: NAME.screens.txt (tools/screen_text.py, splash screens removed by
tools/strip_splash.py); NAME.PRN, the report printed to disk, when the list
was reached. NAME = TEST600_TA_OA_TB_OB_SETTINGS.

  TA   OA  TB    OB   outcome
  150  40  590   40   list, A 41, B 11, 54 pairs; .PRN
  150  40  600   40   list, A 41, B 1, 11 pairs; .PRN ("End: Target + 0")
  150  40  601   40   "Must enter value for this field" (target B)
  150  40  650   40   "Must enter value for this field" (target B)
  150  100 1000  100  "Must enter value for this field" (target B), built-in
  150  100 600   100  "No Matching Primers Found", built-in
  579  40  600   40   list, A 41, B 1, 29 pairs; .PRN
  580  40  600   40   list, A 41, B 1, 30 pairs; .PRN
  581  40  600   40   "Primer too long to fit within molecule." (target A)
  590  40  600   40   "Primer too long to fit within molecule." (target A)
  600  40  600   40   "Primer too long to fit within molecule." (target A)
  601  40  600   40   "Must enter value for this field" (target A)
  620  40  600   40   "Must enter value for this field" (target A)
  650  40  700   40   "Must enter value for this field" (target A)
  0    40  600   40   "Must enter value for this field" (target A)
  576  40  600   40   primer A length 25 (typed "2 5 enter" in its Length
                      field): "Primer too long to fit within molecule."
  5    4   10    40   "Primer too long to fit within molecule."
  5    4   19    40   "Primer too long to fit within molecule."
  5    4   20    40   list, A 5, B 41, 132 pairs
  5    4   21    40   list, A 5, B 41, 127 pairs
The runs from 0 on stop after the form (keys up to "space , ,"); 576 and 0 are
in KEYS as: "... enter 5 7 6 enter 4 0 enter enter 2 5 enter 6 0 0 enter 4 0
enter enter enter enter enter , , , , space , ," and "... enter 0 enter 4 0
enter enter enter 6 0 0 enter 4 0 enter ...", 5/10 (and 19, 20, 21): "... enter 5 enter 4
enter enter enter 1 0 enter 4 0 enter ...".

The four .PRN files equal primerpicker's reports (find-pair --offsets 40,40
--report, with --config PD2_RELAXED.CFG), byte for byte.
Runs that lost keys or froze were rerun; target A 582 froze three times and
is not recorded.
Cursor positions after each box: ../PROBES_20261021.md, section 1.
