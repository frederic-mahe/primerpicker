#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "primerpicker find-pair" when no pair is found.
## Expected values come from PRIMER.exe
## (tests/legacy/search/find_pair_no_match.tsv).
## Usage: test_find_pair_no_match.sh /path/to/primerpicker /path/to/tests/legacy

## Print a header
SCRIPT_NAME="Test find-pair, no match"
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
TLONG="${2}/enzymes/TLONG.GB"

## targets 500 and 1100: primers on both sides, but no pair matches
DESCRIPTION="find-pair: no matching pair is reported"
"${PRIMERPICKER}" find-pair "${TLONG}" 500 1100 | \
    grep -qxF "No Matching Primers Found" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="find-pair: no list header without pairs"
"${PRIMERPICKER}" find-pair "${TLONG}" 500 1100 | \
    grep -q "Position" && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

## targets 120 and 300: no possible primer B
DESCRIPTION="find-pair: no possible primer B is reported as no match"
"${PRIMERPICKER}" find-pair "${TLONG}" 120 300 | \
    grep -qxF "No Matching Primers Found" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit 0
