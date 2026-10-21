Uniqueness and False Priming of complement-strand primers (recorded 2026-10-21)
===============================================================================

Observations: ../PROBES_20261021.md, section 3.

Input (not a recording): TCOMP.GB, written with
  python3 -I tools/plant_sites.py TCOMP 2110 0.5,0.5,0.5,0.5 250 681:12:200 TTCACCGGTTAA@450 689:10:900 CGGTTAAGGGG@300 > TCOMP.GB
1000 bases; the complement-strand primer at 700 is ACAATTCACCGGTTAAGGGG
(bases 681 to 700); copies of parts of it: top strand at 200 (681:12) and
900 (689:10), reverse complements at 300 (CGGTTAAGGGG, with its 3' end) and
450 (TTCACCGGTTAA).

Recorded with (P = position, K = u or f):
  PRIMER_EXE=... PAUSE=3 tools/probe_keys.sh WORKDIR NAME "space enter enter down down
    enter , , enter P P P enter y enter enter enter s , , right right enter , , K , ,
    right , , enter , , down down enter , , r e p o r t enter , , , ," MOLECULE
(the keys of tools/homology_screen.sh for a complement-strand primer, then the
Print menu of the screen, Print to Disk, file REPORT).

Files (NAME = MOLECULE_POSITIONc_K): NAME.screens.txt (splash screens removed
by tools/strip_splash.py), NAME.PRN (the report):
  TREP2_400c_u, TREP2_400c_f   ../homology/TREP2.GB, TGTATCCAGGTAATGCAAGG
  TCOMP_700c_u, TCOMP_700c_f   TCOMP.GB, ACAATTCACCGGTTAAGGGG
  TCOMP_708c_u                 TCOMP.GB, CTAAACCAACAATTCACCGG
  TCSTART_20c_u, TCSTART_28c_u TCSTART.GB (an input: plant_sites.py TCSTART 2111
                               0.5,0.5,0.5,0.5 250 1:12:500 9:12:700),
                               TTAGAGAACTAGGCTTAGTG and the primer at 28
  TREP2_405c_u                 screens only (the report's file name was
                               lost twice; the screen map is recorded)
