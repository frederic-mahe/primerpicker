Mutagenize: a primer ending with a run of 3 (recorded 2026-10-21)
=================================================================

Observations: ../PROBES_20261021_picking.md, "Mutagenize and the run at
the 3' end".

Modify -> Mutagenize on the primer at 450 of ../search/TRUN.GB
(GCAGTCAGGCTCATCGCAAA, top strand, 20 bases), which meets every
criterion of the Primer Analysis screen but ends with AAA; no edit.
Recorded with
  PRIMER_EXE=... PAUSE=3 tools/mutagenize.sh WORKDIR ../search/TRUN.GB 450 ","
File: TRUN_450.screens.txt (splash screens removed by
tools/strip_splash.py). The ORIGINAL and MUTATED lines read
"%GC: 55    Tm°C: 75    Hairpin: OK  Dimers: OK  Runs: !".
