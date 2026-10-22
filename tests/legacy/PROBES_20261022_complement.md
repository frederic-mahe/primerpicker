# Probes of PRIMER.exe, 2026-10-22: complement-strand primers

Two open questions of the proposal (§6.4 and §6.3f): Mutagenize on a
complement-strand primer, and Match Active Primer with a
complement-strand primer A. Recorded with DOSBox-X, headless, `PAUSE`
3.5, 2 or 3 emulators. Recordings and keys: `complement/`
(`README.txt`). Nothing in the README, the tools or the earlier
recordings was changed.

Predictions: the four Match lists were computed before comparison with
primerpicker's `legacy_match_active_primer` given primer A's
complement-strand sequence (a temporary CLI change, not kept): the
lists of `match_c400`, `match_c460_all` and `match_c400_compB` (1, 70
and 7 rows) agree row by row; `match_c400_all` is still to compare. The
Mutagenize molecule line was predicted as the reverse complement of the
molecule (below).


## Draft for tests/legacy/README.md

## Complement-strand primers (recorded 2026-10-22)

`complement/` (`README.txt`): primers A made from TEST600 on the
complement strand at 400(C) and 460(C).

**Mutagenize** accepts them (PRIMER.exe has no restriction to the top
strand). The header line reads `Pos:   400   COMPLEMENTARY`; the box
shows the primer's own strand, 5' to 3': the molecule line is the
reverse complement of the molecule, from 20 bases before the primer's
5' end to 39 after it (`TCATCTGGCG...GCTAAC` for 400(C): positions 181
to 240 of the reverse complement of the 600 bases, the primer's 5' end
being 600 - 400 + 1 = 201), the primer under it, and the reading frames
as for a top-strand primer, on that strand. Typing over bases works as
on the top strand (`mut_c400_edit`: `a` on the third base and `g` on the
eighth: MUTATED %GC 45, Tm 69, Hairpin OK); S shows "Primer Saved"; F2
opens the codon table (`mut_c460`).

**Match Active Primer** with primer A on the complement strand: the
Define Matching Primer form offers "Complement strand (Y/N): N"
(primers B on the top strand; for a top-strand primer A it offers Y);
answered Y, the primers B are on the complement strand too, and the
report reads `Target bp pos:   450 (Complementary strand)`. The list
follows the same rules as for a top-strand primer A, primer A being the
complement-strand sequence: `match_c400` (built-in settings, target
100: 1 row, 73), `match_c400_all` and `match_c460_all`
(`PD2_MATCH_ALL.CFG`), `match_c400_compB` (target 450, primers B on the
complement strand: 7 rows). The report's Primer A line shows the
position without `(C)` (`   400`).
