Find Pair: a key typed during the search (recorded 2026-10-21)
==============================================================

Observations: ../PROBES_20261021_picking.md, "Find Pair's 'Press key to
interrupt search'".

Recorded with tools/find_pair_limits.sh (PAUSE=4), EXTRA being the keys
typed after the form (',' = one more pause, then a space):
  TEST600.GB_300_290_310_290
      PRIMER_EXE=... PAUSE=4 EXTRA=", , , space , , , , , ," SECONDS_MAX=400 \
      tools/find_pair_limits.sh WORKDIR ../pairs/TEST600.GB 300 290 310 290 400
  TEST600.GB_150_100_400_100_PD2_CRITERIA
      PRIMER_EXE=... PAUSE=4 CFG=../config/PD2_CRITERIA.CFG EXTRA=", , space , , , , , , ," \
      SECONDS_MAX=900 tools/find_pair_limits.sh WORKDIR ../pairs/TEST600.GB 150 100 400 100 900
Files: NAME.screens.txt (splash screens removed by tools/strip_splash.py).
In both, the search ends before the space ("<Press key to display
possible pairs>"), which then shows the list.
Not kept: two runs on a 30,000-base random molecule that lost keys.
