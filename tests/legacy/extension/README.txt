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

Later runs (2026-10-22), tools/probe_keys.sh, PAUSE=3.5, SECONDS_MAX=600,
X as above:
  NAME        KEYS                                                    FILES
  ext_unique  X "i , , , u , , , , , ,"                                ../pairs/TEST600.GB
              (Inspect on the picking screen, then U: Uniqueness)
  ext_pair    X "s , , enter , , enter 4 0 0 enter y enter enter enter , , s , ,
              right right down enter 1 enter 2 enter , , a , , , ,"    ../pairs/TEST600.GB
              (primer B at 400(C), Analyze -> Primer Pair 1 2, A: AnalProd)
  ext_multi   "space down down down down down down enter down enter , enter down down enter ,
              right right down down enter 1 2 enter enter , , , , h , , , d , , , ,"
                                                                       P1.PD2 P2.PD2
              (Retrieve from disk P1.PD2 and P2.PD2, Analyze -> Multi Primer 1 2, H, D)
P1.PD2: tools/make_pd2.py GTAAAACGACGGCCAGT GAATTC (= ../primers/M13EXT.PD2);
P2.PD2: tools/make_pd2.py AACAGCTATGACCATG (inputs).
Not kept: ext_false (False Priming; froze before the picking screen).
