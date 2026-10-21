The picking screen: From Molecule and Modify -> Position (recorded 2026-10-21)
===========================================================================

Observations: ../PROBES_20261021_picking.md.

Recorded with
  PRIMER_EXE=... PAUSE=3 tools/probe_picking.sh WORKDIR NAME "$BASE KEYS" MOLECULE
  BASE="space enter enter down down enter , , enter"
(Primer/Create -> From Molecule -> Load File -> the molecule -> Bp
position), then KEYS: the position, Enter on the strand, length and
extension fields (or "y" for the complement strand, which moves on by
itself; or a typed length), then the keys of the picking screen.
probe_picking.sh loads CTRLKEYS.COM (tools/ctrlkeys.S) first: f5 to
f10 stand for Ctrl+Home, Ctrl+End, Ctrl+PgUp, Ctrl+PgDn, Ctrl+Left and
Ctrl+Right. kp_plus and kp_minus are + and -.

Molecules: ../pairs/TEST600.GB, ../search/TRUN.GB, and TWIDTH.GB
(written for these probes, an input: 300 bases, an AT-rich stretch
with one G at 50, a GC-rich stretch at 101 to 200, a mixed one;
python3 -I tools/make_width_molecule.py > TWIDTH.GB).

  NAME             MOLECULE     KEYS
  explore_51       TEST600.GB   "5 1 enter enter enter enter , right , right , right , left , kp_plus , kp_minus , down , up , , ,"
  find_51          TEST600.GB   "5 1 enter enter enter enter , f10 , f10 , f9 , f9 , f9 , , ,"
  pages_51         TEST600.GB   "5 1 enter enter enter enter , pagedown , pageup , home , end , , up , , ,"
  start_2          TEST600.GB   "2 enter enter enter enter , left , left , up , right , , ,"
  end_575          TEST600.GB   "5 7 5 enter enter enter enter , right , right , down , kp_minus , right , kp_plus , kp_plus , , ,"
  end_579          TEST600.GB   "5 7 9 enter enter enter enter , right , right , right , kp_plus , f10 , , ,"
  none_2           TEST600.GB   "2 enter enter enter enter , f9 , , , f10 , , ,"
  save_52          TEST600.GB   "5 1 enter enter enter enter , right , s , , ,"
  ctrl_keys_51     TEST600.GB   "5 1 enter enter enter enter , f5 , f6 , f7 , f8 , , ,"
  find_end         TEST600.GB   "5 3 0 enter enter enter enter , f10 , , , f10 , , , f10 , , ,"
  len_plus         TEST600.GB   "3 6 enter enter enter enter , kp_plus , kp_plus , kp_plus , kp_minus , f10 , , ,"
  comp_400         TEST600.GB   "4 0 0 enter y enter enter , , , right , left , down , f10 , , f9 , , s , , ,"
  esc_52           TEST600.GB   "5 1 enter enter enter enter , right , esc , , , ,"
  modify_position  TEST600.GB   "5 1 enter enter enter enter , s , , right down enter , , , right , f10 , , s , , ,"
  len_99           TEST600.GB   "5 1 enter enter 9 9 enter enter , , kp_plus , , ,"
  w_at_41          TWIDTH.GB    "4 1 enter enter enter enter , , right , , ,"
  w_short_41       TWIDTH.GB    "4 1 enter enter 5 enter enter , , kp_minus , kp_minus , kp_minus , kp_minus , kp_minus , , ,"
  w_gc_120         TWIDTH.GB    "1 2 0 enter enter 4 0 enter enter , , kp_plus , kp_plus , , ,"
  trun_450         TRUN.GB      "4 5 0 enter enter enter enter , , , right , , ,"
  trun_find_40     TRUN.GB      "4 0 enter enter enter enter , f10 , , , f10 , , ,"
  trun_find_440    TRUN.GB      "4 4 0 enter enter enter enter , f10 , , , f10 , , , f9 , , ,"

Files: NAME.screens.txt (tools/screen_text.py, splash screens removed
by tools/strip_splash.py); picking.tsv, every settled state of the
picking screen of these runs (tools/picking_table.py, from the
colours.txt that probe_picking.sh writes, not kept: about 9 MB a run):
run, molecule, frame, strand, primer length, position, MEETS CRITERIA
SELECTED, the value line, the window's background colour. States of
len_99 are not in the table (its window is drawn on 71 bases only).
Not kept: earlier runs that lost keys or froze (esc_52, comp_400, a
first find_51 whose Load File list held CTRLKEYS.COM), and a first
exploration with another tool.
