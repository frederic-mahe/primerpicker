Multi Primer Analysis: the map page (M) (recorded 2026-10-20)
=============================================================

Primers made from a molecule in the same session (Primer/Create -> from
Molecule; the second and third reuse the loaded molecule), then Analyze
-> Multi Primer and M. Recorded with
  PRIMER_EXE=... tools/probe_keys.sh WORKDIR NAME "KEYS" MOLECULE
Each file holds the distinct screens (tools/screen_text.py) from the
Multi Primer Analysis screen on, the map page last.

TEST600_51_300_m.screens.txt: ../pairs/TEST600.GB (600 bases), top
  strand 51 and 300; fields A=1, B=2, C and D empty. PAUSE=2.5,
  KEYS "space enter enter down down enter , , enter 5 1 enter enter
  enter enter s , enter , , enter 3 0 0 enter enter enter enter s ,
  right right down down enter 1 2 enter enter , , m , , pagedown , ,"
TEST600_51_120_400c_m.screens.txt: TEST600.GB, top strand 51 and 120,
  complement strand 400; A=1, B=2, C=3, D empty. PAUSE=2,
  KEYS "space enter enter down down enter , , enter 5 1 enter enter
  enter enter s , enter , enter 1 2 0 enter enter enter enter s ,
  enter , , enter 4 0 0 enter y enter enter enter s , right right down
  down enter 1 2 3 enter , m"
TMAP6600_3300_6500c_m.screens.txt: ../homology/map/TMAP6600.GB (6600
  bases), top strand 3300, complement strand 6500; A=1, B=2. PAUSE=2.5,
  KEYS "space enter enter down down enter , , enter 3 3 0 0 enter enter
  enter enter s , enter , , enter 6 5 0 0 enter y enter enter enter s ,
  right right down down enter 1 2 enter enter , m"

Observations
------------
Multi Primer screen, primers from a molecule: "Primer A (bp    51 ):",
"Primer B (bp  6500C):" (position in 5 columns, then "C" for the
complement strand or a space).
Map page ("* * *   MULTI-PRIMER MOLECULE MAP   * * *"):
- "Source molecule:  NAME" from column 4, then the size starting in
  column 56 (left-aligned, not right-aligned as in the Match report),
  then " bps": TEST600, 27 spaces, "600 bps"; TMAP6600, 26 spaces,
  "6600 bps" (two recordings only);
- "Molecule map with primer locations and directions marked:", then a
  line of '═' from screen column 4 (0-based, the box border being
  column 0) with L / s characters, s = L / 73 + 1 (integer division,
  the scale of the Uniqueness report map, README.md): 66 characters
  for 600 bases (s = 9), 72 for 6600 (s = 91);
- one line per filled field, in the order A to D, under the map: the
  letter in column 4 + P / s (P the primer's position: 51 -> 9, 120 ->
  17, 300 -> 37, 400 -> 48 with s = 9; 3300 -> 40, 6500 -> 75 with s =
  91), followed by '>' (top strand: "A>") or preceded by '<'
  (complement strand: "<C"). (P - 1) / s gives the same columns for
  these six primers;
- "Bp Pos." then "Primer A:" to "Primer D:" with the position right-
  aligned in 8 columns after the colon ("Primer A:      51", "Primer
  A:    3300"), " (C)" for the complement strand, "  none" for an empty
  field;
- bottom line "Press ESCape to return to prior screen"; pagedown
  changes nothing.
