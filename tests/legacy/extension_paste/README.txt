Pasting an enzyme site into a 5' extension (recorded 2026-10-20)
================================================================

Primer made from ../pairs/TEST600.GB at position 51 (top strand, length
20, CTCAAAGGCTCAGACATCAG), built-in settings. Recorded with
  PRIMER_EXE=... tools/probe_keys.sh WORKDIR NAME "KEYS" ../pairs/TEST600.GB
(PAUSE as given); the files hold the distinct screens
(tools/screen_text.py) from the first relevant screen on.

create_form_f1_ecori.screens.txt (PAUSE=3): F1 in the "5' extension to
  be added" field of the Create Primer form.
  KEYS "space enter enter down down enter , , enter 5 1 enter enter
  enter f1 , , e down down down down enter , , enter , , i , , ,"
modify_extension_f1_ecori.screens.txt (PAUSE=2.5): Modify ->
  5'-extension (Change 5' Extension form), F1.
  KEYS "space enter enter down down enter , , enter 5 1 enter enter
  enter enter s , right down down down enter , , f1 , , e down down
  down down enter , , enter , , , ,"
modify_extension_f1_ecori_then_keys.screens.txt (PAUSE=3): the same,
  then "right enter , , , h , , d , ," (meant for Analyze -> Active
  Primer; it opened the Modify Primer Position screen instead, then
  Inspect values are shown through its own keys; kept for the values).
direct_entry_modify_extension.screens.txt (PAUSE=2.5): a Direct Entry
  primer (CAGTCAGACTGACAGTCAGT), then Modify -> 5'-extension.
  KEYS "space down enter enter c a g t c a g a c t g a c a g t c a g t
  enter s , right down down down enter , , f1 , , e , enter , , , ,"

Observations
------------
- F1 opens the same enzyme list as in Mutagenize (name and site, "Use
  Initial letter as fast key    ESC = Exit"); "e" jumps to Eco31, four
  downs to EcoRI; Enter closes the list and writes the site, GAATTC,
  into the extension field ("5' extension to be added:    GAATTC";
  "5' extension:  GAATTC" in the Change 5' Extension form).
- The length (incl 5' extension) stays 20: the primer keeps its 5' end
  (position 51) and loses bases at its 3' end. Primer List after
  Modify -> 5'-extension: "1   TEST600       !     51     20
  GAATTC-CTCAAAGGCTCAGA". Modify Primer Position screen: "Primer
  length: 20   Molecule: TEST600        5'extension: GAATTC", "Pos: 51
  %GC: 45!   Tm°C: 47!   Hairpin: OK  Dimers: OK  Runs: !".
- Inspect (create_form_f1_ecori, key i): "Sequence:
  5'-GAATTCCTCAAAGGCTCAGA-3'", "Length:  20      5' Ext: GAATTC",
  "Composition:  A:7   C:5   G:4   T:4", "From molecule:   TEST600
  bp pos:  51", GC 45 NO, Tm 47 NO, No Hairpins -- YES, 3' dimers 2
  YES, dimers 6 YES, runs 3 NO, 3' GC runs 0 YES. These are the values
  of primerpicker's `analyze --extension GAATTC CTCAAAGGCTCAGA` (Tm of
  the 14-base primer without the extension, the rest on the 20 bases).
- A Direct Entry primer cannot be changed there: Modify -> 5'-extension
  shows a box "Cannot Modify Directly Entered Primers" / "Press any key
  to continue..." (the Direct Entry form has no F1 key).
- A second F1 adds the next site where the cursor is, after the first
  (modify_extension_f1_ecori_f1_bamhi.screens.txt, PAUSE=3, KEYS as
  modify_extension_f1_ecori up to the EcoRI "enter , ,", then "f1 , ,
  b down enter , , enter , , ,"): "5' extension:  GAATTCGGATCC"; Primer
  List "GAATTCGGATCC-CTCAAAGG" (length still 20: 12 + 8). The cursor
  was at the end of GAATTC: whether a site pasted in the middle of an
  extension is inserted or typed over is not recorded.
