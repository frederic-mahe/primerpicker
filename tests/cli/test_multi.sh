#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "primerpicker multi". Expected values come from
## PRIMER.exe screens (tests/legacy/multi/multi_primer_balance.tsv and
## multi_primer_dimers.tsv).
## Usage: test_multi.sh /path/to/primerpicker

## Print a header
SCRIPT_NAME="Test multi"
LINE=$(printf '%.0s-' {1..76})
printf "# %s %s\n" "${LINE:${#SCRIPT_NAME}}" "${SCRIPT_NAME}"

## Declare a color code for test results
RED="\033[1;31m"
GREEN="\033[1;32m"
NO_COLOR="\033[0m"

failure () {
    printf "%bFAIL%b: %s\n" "${RED}" "${NO_COLOR}" "${1}"
    exit 1
}

success () {
    printf "%bPASS%b: %s\n" "${GREEN}" "${NO_COLOR}" "${1}"
}

PRIMERPICKER="${1}"

DESCRIPTION="multi with one primer is an error"
"${PRIMERPICKER}" multi GTAAAACGACGGCCAGT > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="multi with five primers is an error"
"${PRIMERPICKER}" multi AGGAATGGTGATATT AGGAATGGTGATATT AGGAATGGTGATATT \
                  AGGAATGGTGATATT AGGAATGGTGATATT > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="multi with one primer and empty fields is an error"
"${PRIMERPICKER}" multi GTAAAACGACGGCCAGT - - > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="multi rejects letters other than A, C, G and T"
"${PRIMERPICKER}" multi GTAAAACGACGGCCAGT AACAGCTANGACCATG > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

## selection 1234 of multi_primer_*.tsv
MULTI=$("${PRIMERPICKER}" multi CTCAAAGGCTCAGACATCAG GATCCATGCAATGCTAGCTA \
                                ACGTTGCAGGTCAAGCTTGC TTGGCCAAGTCGATCGGATC)

check_line () {
    local DESCRIPTION="multi: ${1}"
    grep -qxE "${2}" <<< "${MULTI}" && \
        success "${DESCRIPTION}" || \
            failure "${DESCRIPTION}"
}

check_line "primer D" "Primer D: +5'-TTGGCCAAGTCGATCGGATC-3'"
check_line "length" "Length +20 +20 +20 +20"
check_line "GC %" "% GC +50 +45! +55 +55"
check_line "Tm" "Tm °C +66 +67 +74 +76"
check_line "3' G/C runs" "Runs of G/C at 3' end +1 +0 +2 +1"
check_line "runs" "Runs of bases \(any\) +3 +2 +2 +2"
check_line "hairpins" "Hairpins +OK +OK +! +OK"
check_line "3' dimers, row B" "  B +\. +6! +4! +1"
check_line "3' dimers, row D" "  D +\. +\. +\. +4!"
check_line "any dimers, row C" "  C +\. +\. +8! +4"

## selection 1_2: field B left empty
DESCRIPTION="multi: an empty field leaves its column empty"
"${PRIMERPICKER}" multi CTCGATCCGTTCCTAAT - CGGGGCTCATTGCACTC | \
    grep -qxE "Primer B:" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit 0
