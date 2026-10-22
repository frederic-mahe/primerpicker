Complement-strand primers: Mutagenize and Match Active Primer (recorded 2026-10-22)
====================================================================================

Observations: ../PROBES_20261022_complement.md.

Primer A made from ../pairs/TEST600.GB on the complement strand (20
bases, built-in settings unless PD2.CFG), recorded with
  PRIMER_EXE=... PAUSE=3.5 [SECONDS_MAX=600] tools/probe_keys.sh WORKDIR NAME "KEYS" ../pairs/TEST600.GB [PD2.CFG]
with
  A(POS) = "space enter enter down down enter , , enter POS enter y enter enter , , s , ,"
           (From Molecule -> Load File -> TEST600.GB, Bp position POS,
           Complement strand y, picking screen, Save+Exit)
  M(MOL, POS, TARGET, STRAND) = "space enter enter MOL enter , , enter POS enter y enter enter , , s , ,
           down down enter , TARGET enter enter enter STRAND enter enter enter enter , , ,
           p , down down enter , r e p o r t enter , , ,"
           (then Match Active Primer, TARGET, the Complement strand field
           answered STRAND, Print to Disk as REPORT)

  NAME              SECONDS_MAX  KEYS
  mut_c400          -            A("4 0 0") " right down down down down down enter , , , , ,"
  mut_c400_edit     600          A("4 0 0") " right down down down down down enter , , , right right a right right right right g , , s , , , ,"
  mut_c460          600          A("4 6 0") " right down down down down down enter , , , f2 , , , ,"
  match_c400        -            A("4 0 0") " down down enter , 1 0 0 enter enter enter enter enter enter enter enter enter , , , p down down enter r e p o r t enter , , ,"
  match_c400_all    600          M("down down down", "4 0 0", "1 0 0", enter), PD2.CFG = ../match_variations/PD2_MATCH_ALL.CFG
  match_c460_all    600          M("down down down", "4 6 0", "3 0 0", enter), PD2.CFG = ../match_variations/PD2_MATCH_ALL.CFG
  match_c400_compB  600          M("down down", "4 0 0", "4 5 0", y)

Files: NAME.screens.txt (tools/screen_text.py, splash screens removed by
tools/strip_splash.py); NAME.PRN: the REPORT file as written (DOS line
ends, code page 437).
Not kept: runs that froze early (3 emulators at once).
