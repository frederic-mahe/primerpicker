Mutagenize: F1 (enzyme list) and S (save) (recorded 2026-10-20)
=================================================================

Modify -> Mutagenize on the primer made from ../pairs/TEST600.GB at
position 51 (CTCAAAGGCTCAGACATCAG, top strand, 20 bases), recorded with
  PRIMER_EXE=... tools/mutagenize.sh WORKDIR ../pairs/TEST600.GB 51 "KEYS"
(PAUSE 3, the tool's default). Files:
- NAME.screens.txt: the distinct screens (tools/screen_text.py) from
  the first Mutagenize Primer screen on. screen_text.py prints each
  distinct screen once: a screen identical to an earlier one (for
  instance the Mutagenize screen after the enzyme list is closed) is
  not printed again, so the last screen of the file is not always the
  last state;
- NAME.final.png: the last frame of the recording (the last state);
- TEST600_51_f1_pagedown_pagedown_end.list_states.png: the enzyme list
  after F1, after the first pagedown, after the second (crops of three
  frames, left to right; the highlight is the grey bar).

KEYS of each NAME ("GAATTC" = g right a right a right t right t right c,
which types GAATTC over the first six bases):
  TEST600_51_f1_pagedown_pagedown_end    "f1 , , pagedown , pagedown , end ,"
  TEST600_51_GAATTC_f1_pagedown_pagedown "GAATTC , f1 , , pagedown , pagedown ,"
  TEST600_51_f1_enter                    "f1 , enter , , ,"
  TEST600_51_right3_f1_e_enter           "right right right f1 , e , enter , , ,"
  TEST600_51_GAATTC_f1_e_down4_enter     "GAATTC , f1 , e down down down down enter , , ,"
  TEST600_51_GAATTC_s_enter              "GAATTC , s , , , enter , , ,"
  TEST600_51_GAATTC_s_enter_esc          "GAATTC , s , enter , esc , , ,"

Observations
------------
F1 opens a box over the right part of the screen: the enzyme list,
name and recognition site, 14 rows per page, starting "AatII GACGTC",
"AflII CTTAAG", "AluI AGCT", ... with a scroll bar; the bottom line
reads "Use Initial letter as fast key    ESC = Exit". The list is the
same before and after an edit that creates an EcoRI site (no marking of
sites in the primer). Keys: pagedown moves the highlight to the last
row of the page (BbeI), a second pagedown shows the next page, which
starts with that row (BbeI to DraI, DraI highlighted); end changes nothing; a letter jumps to the
first enzyme with that initial ("e": Eco31, then four downs: EcoRI).
Enter closes the list. In the three runs (AatII with the edit
highlight on base 1; Eco31 on base 4; EcoRI after the GAATTC edit,
highlight on base 6) nothing on the screen changed: the primer, the
translations and the MUTATED line are those from before F1
(NAME.final.png).
S (after the GAATTC edit; MUTATED line: GC 50, Tm 68, Hairpin OK,
Dimers OK, Runs OK) shows a box "Primer Saved" / "Press any key to
continue..." over the Mutagenize screen; a key closes it and the
Mutagenize screen stays. After Esc, the Primer List holds the original
primer (1, TEST600, Crit "!", Pos 51) and a new entry 2: "TEST600, OK,
51, 20, GAATTCGGCTCAGACATCAG", highlighted (active), with a mark drawn
in the column after its number (a glyph screen_text.py does not decode;
it looks like the code page 437 character 0x7F, a small house:
TEST600_51_GAATTC_s_enter_esc.final.png).
Not known: what the enzyme chosen with Enter is used for (if anything);
what the mark of entry 2 means.
