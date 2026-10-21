Hairpins with tandem G.T pairs (recorded 2026-10-21)
====================================================

Observations: ../PROBES_20261021.md, section 6 (primer B 1867 of
search/TRUN.GB, ATCTGTGGGACCGCCTGCTT, whose screens are in
../find_pair_trun/INSPECT_SCREENS_ATCTGTGGGACCGCCTGCTT.screens.txt).

Two primers designed with a 6-pair stem holding two adjacent G.T pairs,
a loop of AAAA, and flanks, typed with Direct Entry (built-in settings):
  PRIMER_EXE=... PAUSE=3 tools/inspect_primer.sh WORKDIR PRIMER
Files: PRIMER.screens.txt (Primer Analysis, Hairpins and Dimers screens;
splash screens removed by tools/strip_splash.py).

  ACCGTGGAAAACCGTGGTTA  stem CCGTGG / GGTGCC: C:G C:G G.T T.G G:C G:C
                        (tandem G.T then T.G, 5' side first)
                        PRIMER.exe: Hairpin Energy = -5.5 kcals,
                                     5'-ACCGTGGAA
                                        :::::::  )
                                   3'-ATTGGTGCCAA
                        (primerpicker analyze: -5.5)
  ACCTGGCAAAAGCTGGGTTA  stem CCTGGC / GGGTCG: C:G C:G T.G G.T G:C C:G
                        (tandem T.G then G.T, as in ATCTGTGGGACCGCCTGCTT)
                        PRIMER.exe: Hairpin Energy = -6.1 kcals,
                                     5'-ACCTGGCAA
                                        :::::::  )
                                   3'-ATTGGGTCGAA
                        (primerpicker analyze: -0.1, another hairpin)
Both drawings show 7 pairs: the outer A:T pair is part of the stem.
