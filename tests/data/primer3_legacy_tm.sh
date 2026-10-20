#!/bin/sh
# Write the Tm that Primer3's oligotm gives, with the settings that
# reproduce Primer Designer 2.0, for the primers recorded in
# tests/legacy/tm/direct_entry_tm.tsv (primers of A, C, G and T without
# 5' extension, up to 36 bases: oligotm's limit).
#
# usage: sh primer3_legacy_tm.sh > primer3_legacy_tm.tsv
#
# Settings: Breslauer et al. 1986 parameters and Schildkraut salt
# correction (-tp 0 -sc 0), no divalent cations or dNTP (-dv 0 -n 0),
#   DNA:  C = 4e9 exp((10.8 - b) / 1.987) nM with b = 60.6997 cal/K/mol
#         (Primer Designer's constant, tests/legacy/README.md, "Primer
#         Tm"): 0.049611885111237056 nM, so that 10.8 - R ln(C/4) = b
#   salt: 1000 x 10^(0.15 / 16.6) mM = 1021.0244559195794 mM, so that
#         Primer3's 16.6 log10([K+]) - 273.15 equals Primer Designer's
#         -273
set -eu

legacy=$(cd "$(dirname "$0")/../legacy" && pwd)
# oligotm does not print its version: write it in the header by hand
# (the table was made with the Ubuntu package primer3 2.6.1-6build1)
printf '# primer\ttm_celsius (oligotm -tp 0 -sc 0 -dv 0 -n 0 -d 0.049611885111237056 -mv 1021.0244559195794; Primer3 2.6.1)\n'
grep -v '^#' "$legacy/tm/direct_entry_tm.tsv" | \
    awk -F'\t' '$1 == "-" && $2 ~ /^[ACGT]+$/ && length($2) <= 36 {print $2}' | \
    sort -u | \
    while read -r primer; do
        tm=$(oligotm -tp 0 -sc 0 -dv 0 -n 0 -d 0.049611885111237056 -mv 1021.0244559195794 "$primer")
        printf '%s\t%s\n' "$primer" "$tm"
    done
