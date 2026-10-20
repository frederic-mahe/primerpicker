Match Active Primer: variations of the report (recorded 2026-10-20)
======================================================================

Primer/Create -> Match Active Primer (Find Matching Primer list),
printed to disk (P -> Print Control -> Print to Disk, file REPORT).
PRIMER.exe run in DOSBox-X 2026.01.02, headless; the .PRN files are the
REPORT files as written (DOS line ends, code page 437: byte F8 is the
degree sign of "Tm°C").

Inputs
------
PD2_*.CFG: settings files, copied as PD2.CFG next to PRIMER.EXE.
Written with tools/make_pd2_cfg_match.py (built-in defaults, then):
  PD2_MATCH_GC10_TM20.CFG       --gc 10 --tm 20
  PD2_MATCH_100_150.CFG         --start 100 --end 150
  PD2_MATCH_0_10.CFG            --start 0 --end 10
  PD2_MATCH_200_200.CFG         --start 200 --end 200
  PD2_MATCH_ALL.CFG             --gc 100 --tm 99 --criteria 0,100,0,99,-99,99,99,99,99
  PD2_MATCH_ALL_LENGTH25.CFG    the same, --length 25
  PD2_MATCH_ALL_200.CFG         --gc 100 --tm 99 --start 200 --end 200 --criteria (as ALL)
  PD2_MATCH_ALL99_49_49.CFG     --gc 99 --tm 99 --start 49 --end 49 --criteria (as ALL)
  PD2_MATCH_ALL_TOL0_60_60.CFG  --gc 0 --tm 0 --start 60 --end 60 --criteria (as ALL)
  PD2_OFFSETS_11_22_33_44.CFG   --offsets-raw 11,22,33,44 (int16 values 11, 22,
                                33, 44 at file offsets 14, 16, 18, 20)
--start writes file offset 14 and --end offset 16 (the order of
README.md's table); the program reads them the other way round (see
"Settings file" below), so the names above give the bytes, not what the
program shows: PD2_MATCH_100_150.CFG prints "Start: Target - 150 bases;
End: Target + 100 bases", PD2_MATCH_0_10.CFG "Start: Target - 10 bases;
End: Target + 0 bases".
Molecules: ../homology/map/TMAP50K.GB (50000 bases),
../homology/map/TMAP6600.GB (6600), ../search/TRUN.GB (2000).

Recordings
----------
NAME.PRN, NAME = MOLECULE_POSITION_TARGET[_CFG]: primer A (top strand,
20 bases, at POSITION) made from the molecule, then Match Active Primer
with TARGET and the form's other fields left as they are (defaults of
the CFG, built-in without one). Recorded with
  PRIMER_EXE=... tools/match_report.sh WORKDIR MOLECULE POSITION TARGET [CFG]
(PAUSE=2.5 or 3), except:
  TMAP6600_1234_1700.PRN: tools/print_report.sh WORKDIR match
      ../homology/map/TMAP6600.GB 1234 1700 (PAUSE=2.5);
  TRUN_50_600_top_strand.PRN: primer B on the top strand ("Complement
      strand (Y/N)" answered n), recorded with tools/probe_keys.sh,
      keys "space enter enter down down enter , , enter 5 0 enter enter
      enter enter s , down down enter 6 0 0 enter enter enter n enter
      enter enter enter , , p , down down enter , r e p o r t enter , ,"
      and ../search/TRUN.GB.
*.screens.txt: screens (tools/screen_text.py), from the frame before the
search on, for runs that print nothing:
  TRUN_50_30_no_match.screens.txt (print_report.sh match
      ../search/TRUN.GB 50 30, built-in settings),
  TMAP6600_123_400_no_match.screens.txt (print_report.sh match
      ../homology/map/TMAP6600.GB 123 400; primer A has GC 65 %),
  TRUN_50_1000_PD2_MATCH_ALL_200_too_many.screens.txt
      (match_report.sh ../search/TRUN.GB 50 1000 PD2_MATCH_ALL_200.CFG),
  PD2_OFFSETS_11_22_33_44.screens.txt: Modify -> System Defaults ->
      Defaults with PD2_OFFSETS_11_22_33_44.CFG as PD2.CFG
      (tools/probe_keys.sh, keys "space right down down down down down
      down down enter , enter , , esc , ,").
Runs that lost keys (not kept): TMAP50K 950 12000 (keys lost after s),
TMAP6600 300 700 and 5 300 (no matching primer: primer A has GC 35 and
40 %), and a first ALL_200 run.

Observations
------------
Field widths (all 16 reports):
- "Source molecule:  NAME", the size right-aligned to end in column 57,
  then " bps" (line length 61 for 4- and 5-digit sizes).
- "   Target bp pos:" then the target in 6 columns, right-aligned
  ("pos:    30", "pos:   500", "pos:  1700", "pos: 10500"), then
  " (Complementary strand)" for a complement-strand primer B, nothing
  for a top-strand one ("   Target bp pos:   600").
- "   Start: Target - S bases;   End: Target + E bases": S and E without
  padding (1 to 5 digits seen).
- "   GC Match:  +/- %2d GC %      Tm Match:  +/- %2d C".
- Primer A row: position right-aligned in columns 1-6 ("    50",
  "   937", "  1234", " 10050").
- Primer B rows: position right-aligned in columns 1-8 ("      20",
  "     475", "    1675", "   10528"), "bps from target" in columns
  9-18 ("        28", "       162").
- The report ends after the last row (CR LF); no form feed and no page
  break, even with 99 rows (TRUN_50_1000_PD2_MATCH_ALL99_49_49.PRN,
  123 lines).
Search window:
- complement strand: rows from max(target - start offset, primer
  length) to min(target + end offset, molecule length), position =
  5' end; the report prints the offsets as cut: target 30, length 20:
  "Start: Target - 10" (first row 20); length 25: "Start: Target - 5"
  (first row 25); target 1990 of 2000: "End: Target + 10" (last row
  1992; the primers B at 1993 to 2000 have GC 45 % or less, outside
  primer A's 60 +/- 10).
- target 10 (below the primer length): "Start: Target - 65526 bases"
  (10 - 20 = -10 printed as an unsigned 16-bit number), rows from 20.
- top strand (TRUN_50_600_top_strand.PRN): rows from 552 to 609 in the
  window 550-650; "bps from target" = |position - target|.
Tolerances: GC 100 in the settings file is shown as "10" in the form's
two-character field and applied as 10 (TRUN_*_PD2_MATCH_ALL.PRN: "GC
Match:  +/- 10 GC %", rows within 10 % of primer A); 99 is applied as
99.
No matching primer: a box "No Matching Primers Found" / "Press any key
to continue..." over the Define Matching Primer form; a key returns to
the form; nothing is printed.
More than 100: with every position of a 401-position window matching
(PD2_MATCH_ALL_200.CFG), a box "More than 100 sites found.  Suggest you
redefine matching limits." / "Press any key to continue..." over the
Find Matching Primer list screen (rows drawn behind it); a key returns
to the Define Matching Primer form; nothing is printed. With 121
positions that meet the criteria but 1 within the tolerances
(PD2_MATCH_ALL_TOL0_60_60.CFG) the list (1 row) is shown and printed:
the limit counts rows within the tolerances. 99 rows are listed and
printed (PD2_MATCH_ALL99_49_49.CFG).
Settings file: PD2_OFFSETS_11_22_33_44.CFG shows "Match Active Primer:
Start @ - 22 / End @ + 11" and "Find Pair (normal) Start @ - 44 /
(complement) End @ + 33": offset 14 holds the Match Active Primer end
(+), 16 its start (-), 18 the Find Pair complement end (+), 20 the Find
Pair normal start (-). The Match reports agree for offsets 14 and 16.
