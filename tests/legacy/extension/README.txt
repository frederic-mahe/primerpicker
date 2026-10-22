Primers with a 5' extension (recorded 2026-10-22)
=================================================

Observations: ../PROBES_20261022_extension.md.

Primer/Create -> From Molecule on ../pairs/TEST600.GB at 51, length 20
(extension included), 5' extension GAATTC, recorded with
  PRIMER_EXE=... PAUSE=3.5 SECONDS_MAX=600 tools/probe_picking.sh WORKDIR NAME "X KEYS" ../pairs/TEST600.GB  (ext_picking3)
  PRIMER_EXE=... PAUSE=3.5 SECONDS_MAX=600 tools/probe_keys.sh WORKDIR NAME "X KEYS" ../pairs/TEST600.GB     (ext_length2)
  X = "space enter enter down down enter , , enter 5 1 enter enter enter g a a t t c enter , , ,"

  NAME          KEYS
  ext_picking3  "right , , kp_plus , , kp_minus , , f10 , , , s , , , ,"
                (the picking screen: Right, +, -, Ctrl+Right, S)
  ext_length2   "s , , , right down down enter , , , 2 5 enter , , , , ,"
                (S, then Modify -> Length: 25)

Files: NAME.screens.txt (tools/screen_text.py, splash screens removed
by tools/strip_splash.py); picking_extension.tsv: the settled states of
ext_picking3 (tools/picking_table_extension.py pairs/TEST600.GB
WORKDIR/ext_picking3; its colours.txt is not kept, about 9 MB). Not
kept: earlier runs that froze early or lost keys (ext_picking,
ext_picking2, ext_length; "plus" is not an autotype key: kp_plus).
