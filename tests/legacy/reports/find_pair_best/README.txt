Find Pair reports with a typed "Print best: [N]" (recorded 2026-10-21)
======================================================================

Observations: ../../PROBES_20261021.md, section 2.

Recorded with tools/print_report.sh find_pair and PRINT_KEYS (PAUSE 3, or 3.5
for the 4-digit runs):
  TFP_150_400_bestN (N = 25, 30, 50, 0):
    PRIMER_EXE=... PRINT_KEYS=", , p , , D D enter , , down down enter , , r e p o r t enter , , , ," \
      tools/print_report.sh WORKDIR find_pair search/TFP.GB 150 400
    (D D = the digits of N; "0 enter" for 0)
  TFP_150_400_best9999:
    PRINT_KEYS=", p , , 9 , 9 , 9 , 9 , , down down enter , r e p o r t enter , ,"
  TEST600_150_400_CRITERIA_best100:
    CFG=config/PD2_CRITERIA.CFG PRINT_KEYS=", , p , , 1 0 0 enter , , down down enter , , r e p o r t enter , , , ," \
      tools/print_report.sh WORKDIR find_pair pairs/TEST600.GB 150 400
  TEST600_150_400_CRITERIA_best1000:
    CFG=config/PD2_CRITERIA.CFG PRINT_KEYS=", p , , 1 , 0 , 0 , 0 , , down down enter , r e p o r t enter , ,"
    (a fourth digit completes the field by itself: no Enter; two earlier runs
    with an Enter after the fourth digit chose Abandon Printing and wrote
    nothing)

Files: NAME.screens.txt (splash screens removed by tools/strip_splash.py) and
NAME.PRN (the file written, as it is). best0 wrote no file.
