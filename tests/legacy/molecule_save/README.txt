Molecule-Ops -> Save/Remove (S) and Create New (recorded 2026-10-22)
====================================================================

Observations: ../PROBES_20261022_molecule_save.md.

Recorded with
  PRIMER_EXE=... PAUSE=P SECONDS_MAX=M tools/probe_keys.sh WORKDIR NAME "KEYS" FILE...
(M: the default of probe_keys.sh when not given), with these key prefixes:
  L  = "space enter enter down down enter , , esc , , right right right down down enter , ,"
       (From Molecule -> Load File -> the molecule, Esc, then Molecule-Ops -> Save/Remove)
  V  = "space right right right down enter , enter , down down enter , , ,"
       (Molecule-Ops -> View/Edit -> Load File -> the molecule)
  C1 = "space right right right down enter , down enter , z n enter 1 2 enter enter enter e enter b enter d e enter ,"
       (View/Edit -> Create New: name zn, size 12, DNA, linear, left end e,
       right end b, description de)
  C2 = "space right right right down enter , down enter , z n enter 1 2 enter"
  C3 = C2 " enter enter enter enter enter enter , y , , insert , , g o enter d o enter"
       (name zn, size 12, the other fields left empty, Enter Features: y,
       INS, feature name go, description do)

  NAME            P    M    KEYS                                                     FILE
  save_s          3.5  -    L " s , , , enter , , , ,"                               ../pairs/TEST600.GB
  save_new        3.5  -    L " s , , , z z period g b enter , , , , ,"              ../pairs/TEST600.GB
  save_fa         3.5  -    L " s , , , enter , , y , , , ,"                         ../molecules/TFASTA.FA
  save_circ       3.5  -    L " s , , , z c period g b enter , , , , ,"              ../circular/TCIRC.GB
  reload_header   3.5  -    V " down down down enter , , , , ,"                      written/save_new.ZZ.GB as ZZ.GB
  reload_primer   3.5  -    "space enter enter down down enter , , enter , , , ,"    written/save_new.ZZ.GB as ZZ.GB
  create_open     3.5  -    "space right right right down enter , down enter , , , , , ,"
  create_nothing  3    500  C1 " n , , n , , y , , z m enter , , , ,"
  create_feat4    3    500  C1 " y , , insert , , g o enter d o enter enter 5 enter 9 enter , , esc , , n , , y , , z n enter , , , ,"
  create_seq4     3    500  C1 " n , , y , , , , , , , , esc , , , , , ,"
  create_prot     3    500  C2 " p y enter enter enter , , n , , n , , y , , z p enter , , , ,"
  feat_marker4    3.5  700  C3 " m enter 5 enter , y enter , , esc , , esc , , , n , , , y , , , z f enter , , ,"
  feat_region2    3.5  700  C3 " r enter 5 enter 9 enter , , esc , , , n , , , y , , , z f enter , , ,"
  feat_gback      3.5  700  C3 " enter 9 enter 5 enter , , esc , , , n , , , y , , , z f enter , , ,"
  reload_zn       3    500  V " enter , , , , esc , , down down down enter , , , ,"   written/create_feat4.ZN as ZN
  hand1 .. hand4  3    500  V " enter , , , esc , , down down enter , , , esc , , down enter , , , ,"   inputs/ZH1 .. ZH4
  hand6, hand7    3.5  700  V " enter , , , esc , , down down down enter , , ,"     inputs/ZH6, ZH7
  hand8, hand9,   3.5  600  V " enter , , , , ,"                                     inputs/ZH8, ZH9, ZHA
  handA
  tfeat_view      3.5  700  V " enter , , , , , ,"                                   inputs/TFEAT.GB
  round_zh6       3.5  600  L " s , , , enter , , y , , , ,"                         inputs/ZH6
  round_zh7       3.5  600  L " s , , , enter , , y , , , ,"                         inputs/ZH7
                  (loaded with From Molecule, then saved over itself: written/round_zh6.ZH6,
                  written/round_zh7.ZH7; recorded later on 2026-10-22)

Files:
- NAME.screens.txt: tools/screen_text.py, splash screens removed by
  tools/strip_splash.py;
- written/NAME.FILE: the file PRIMER.exe wrote in run NAME (save_fa
  overwrote TFASTA.FA, after "File already exists!  Overwrite (Y/N)?");
- inputs/: files written for these probes (inputs, not recordings):
  ZH1 to ZHA in the native format, by hand (ZH1 to ZH4 with printf from a
  shell, which dropped the NUL bytes: their feature lines read "gada",
  "gbdb", "gcdc"); TFEAT.GB is ../pairs/TEST600.GB with three FEATURES
  entries.
