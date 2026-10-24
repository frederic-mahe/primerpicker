# Example molecule

`pUC19.gb` is the molecule of Help -> Example Molecule (and of
`primerpicker-gui --demo`), built into the graphical program: the
cloning vector pUC19, 2686 bp, circular, GenBank accession L09137
(version L09137.2), from

> Yanisch-Perron, C., Vieira, J., Messing, J. (1985) Improved M13 phage
> cloning vectors and host strains: nucleotide sequences of the M13mp18
> and pUC19 vectors. Gene 33: 103-119.

It was chosen (2026-10-24) because lab users know it: it holds the M13
forward (-20) and reverse primer sites (bases 379 and 461 of the
record, the reverse primer on the complement strand) and the
polylinker of restriction sites (EcoRI at 396, HindIII at 447).

The file is the record as NCBI sends it, unchanged, downloaded on
2026-10-24 with

```sh
curl -o pUC19.gb 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=L09137&rettype=gbwithparts&retmode=text'
```

SHA-256: `9498daed1c3703dd739c684a3a7bcd6770823aa4a500cf9b57cf59bbdc38f72a`

NCBI places no restrictions on the use or distribution of the GenBank
data (<https://www.ncbi.nlm.nih.gov/home/about/policies/>).
