Molecule-Ops -> View/Edit and Save/Remove (recorded 2026-10-21)
===============================================================

Observations: ../PROBES_20261021_picking.md, "Molecule-Ops -> View/Edit".

Recorded with
  PRIMER_EXE=... PAUSE=P tools/probe_picking.sh WORKDIR NAME "KEYS" ../pairs/TEST600.GB
V standing for "space right right right down enter , enter , down down enter , , ,"
(Molecule-Ops -> View/Edit -> Load File -> TEST600.GB):

  NAME             P    KEYS
  view_open        3    "V , , ,"
  view_down        3    "V down , enter , , , , ,"                  (Sequence)
  seq_edit         3    "V down enter , , enter , , , ,"            (Sequence, then Edit)
  features         3    "V enter , , , , ,"
  annotations      3    "V down down enter , , , , ,"
  header           3    "V down down down enter , , , , ,"
  seq_search       3    "V down enter , , right right , , enter , , , ,"
  find_catcta      3.5  "V down enter , , , right right , , enter , , c a t c t a enter , , , ,"
  find_ttaaat      3.5  "V down enter , , , right right , , enter , , t t a a a t enter , , , ,"
  find_catnta      3.5  "V down enter , , , right right , , enter , , c a t n t a enter , , , ,"
  find_aaat        4    "V down enter , , , right right , , enter , , a a a t enter , , , ,"
  saveremove_open  3    "space right right right down down enter , enter , down down enter , , , , , ,"
                        (Molecule-Ops -> Save/Remove, no molecule loaded)
  saveremove_loaded 3.5 "space enter enter down down enter , , esc , , right right right down down enter , , , , , ,"
                        (TEST600 loaded with From Molecule, then Save/Remove)

Files: NAME.screens.txt (tools/screen_text.py, splash screens removed by
tools/strip_splash.py). Not kept: runs of the Format menu (Double
Strand, Translate, Block Size Toggle): neither arrows nor initial
letters changed the view in five runs.
