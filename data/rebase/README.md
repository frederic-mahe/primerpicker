# REBASE data

The restriction enzymes of PrimerPicker come from REBASE, The
Restriction Enzyme Database (http://rebase.neb.com), version 610
(September 2026), in the formats made for EMBOSS:

| File | Contents |
|---|---|
| `emboss_e.610` | name, recognition sequence, cut positions |
| `emboss_r.610` | organism, isoschizomers, suppliers, references |
| `emboss_s.610` | supplier codes |
| `withrefm.610` | all REBASE enzymes: recognition sequences with cuts, isoschizomers, suppliers, references |

The files are copied unchanged from ftp://ftp.neb.com/pub/rebase/
(downloaded 2026-10-18; SHA-256 below).

REBASE is a free and independent resource. Its documentation
(`REBASE.DOC`) allows distributing REBASE files with software,
provided users are not charged for the REBASE data and the following
reference is cited:

> Roberts, R.J., Vincze, T., Posfai, J., Macelis, D. (2010)
> REBASE--a database for DNA restriction and modification: enzymes,
> genes and genomes. Nucl. Acids Res. 38: D234-D236.

PrimerPicker includes the commercially available enzymes (those with
at least one supplier in `emboss_r.610`). `make_enzyme_table.py`
writes them as a C++ table, `src/core/enzymes/rebase_table.hpp`:

```sh
python3 -I data/rebase/make_enzyme_table.py data/rebase 610 > src/core/enzymes/rebase_table.hpp
```

The other REBASE enzymes with a known recognition sequence (4240 in
version 610, from `withrefm.610`) are accepted as names only: old
names (Asp718 is Asp718I) and enzymes no longer sold. Their table is
`src/core/enzymes/rebase_other_table.hpp`:

```sh
python3 -I data/rebase/make_enzyme_table.py data/rebase 610 others > src/core/enzymes/rebase_other_table.hpp
```

Cuts in `withrefm` are converted to the numbering of the EMBOSS files
(G^AATTC: 1 5; GAAGAC(2/6): 8 12; ^GATC: -1 4); the conversion gives
the `emboss_e.610` values for 578 of the 581 commercial enzymes (the
other three are written with trailing N, as CASTGNN^, whose bottom
cut `withrefm` does not give: such cuts are left unknown).

To update REBASE (once a year, see the maintenance calendar), download
the four files of the new version, run both commands with the new
version number, and review the differences of the generated tables.

SHA-256:

```
3e11a384987203c30554504420cd35c18a0fa83894cd737213e5eb29986c29a3  emboss_e.610
4bfd69408bd2b0e19c195be5032b5c563fd5f7071b112dd95e44ce59f84ececa  emboss_r.610
a34f6cd5b04a2ff6e6b38574eab52c296b23a33a4fd1aa9003c0e38f86ea3df5  emboss_s.610
a12c968e421363de250016006195391e781618bb66d9732d6948d47b0b60dd7f  withrefm.610
```
