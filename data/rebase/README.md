# REBASE data

The restriction enzymes of PrimerPicker come from REBASE, The
Restriction Enzyme Database (http://rebase.neb.com), version 610
(September 2026), in the formats made for EMBOSS:

| File | Contents |
|---|---|
| `emboss_e.610` | name, recognition sequence, cut positions |
| `emboss_r.610` | organism, isoschizomers, suppliers, references |
| `emboss_s.610` | supplier codes |

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

To update REBASE (once a year, see the maintenance calendar), download
the three files of the new version, run the command with the new
version number, and review the differences of the generated table.

SHA-256:

```
3e11a384987203c30554504420cd35c18a0fa83894cd737213e5eb29986c29a3  emboss_e.610
4bfd69408bd2b0e19c195be5032b5c563fd5f7071b112dd95e44ce59f84ececa  emboss_r.610
a34f6cd5b04a2ff6e6b38574eab52c296b23a33a4fd1aa9003c0e38f86ea3df5  emboss_s.610
```
