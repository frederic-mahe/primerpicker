#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "primerpicker match --top-strand": primers B on the
## top strand. Expected values come from PRIMER.exe
## (tests/legacy/match_variations/TRUN_50_600_top_strand.PRN, recorded
## 2026-10-20).
## Usage: test_match_top_strand.sh /path/to/primerpicker /path/to/tests/legacy

## Print a header
SCRIPT_NAME="Test match, top strand"
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
LEGACY="${2}"
TRUN="${LEGACY}/search/TRUN.GB"
RECORDED="${LEGACY}/match_variations/TRUN_50_600_top_strand.PRN"
WORKDIR=$(mktemp -d)
trap 'rm -rf "${WORKDIR}"' EXIT
REPORT="${WORKDIR}/REPORT"

## the report without its date (bash only: no cmp in MSYS2)
body () {
    sed 's|Report Date: [0-9/]*|Report Date: DATE|' "${1}"
}

DESCRIPTION="match --top-strand --report: PRIMER.exe's report, but for the date"
"${PRIMERPICKER}" match --top-strand --report "${REPORT}" "${TRUN}" 50 600 > /dev/null && \
    [[ "$(body "${REPORT}"; printf x)" == "$(body "${RECORDED}"; printf x)" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="match --top-strand: the screen names the top strand"
"${PRIMERPICKER}" match --top-strand "${TRUN}" 50 600 | \
    grep -q "primers B on the top strand" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="match --top-strand: the options in either order"
"${PRIMERPICKER}" match --report "${REPORT}" --top-strand "${TRUN}" 50 600 > /dev/null && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="match: an option given twice is an error"
"${PRIMERPICKER}" match --top-strand --top-strand "${TRUN}" 50 600 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

exit 0
