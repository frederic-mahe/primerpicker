#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "primerpicker find-pair" at the 1000-pair limit:
## the best 1000 of the first 1001 pairs. Expected values come from
## PRIMER.exe (tests/legacy/README.md, "Reports of more than 20 pairs";
## tests/legacy/find_pair_trun/, recorded 2026-10-21).
## Usage: test_find_pair_cap.sh /path/to/primerpicker /path/to/tests/legacy

## Print a header
SCRIPT_NAME="Test find-pair, 1000-pair limit"
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
TEST600="${LEGACY}/pairs/TEST600.GB"
TRUN="${LEGACY}/find_pair_trun"
WORKDIR=$(mktemp -d)
trap 'rm -rf "${WORKDIR}"' EXIT
REPORT="${WORKDIR}/REPORT"

## the report without its date (bash only: no cmp in MSYS2)
body () {
    sed 's|Report Date: [0-9/]*|Report Date: DATE|' "${1}"
}

## weights 1 0 1 3 1, targets 260 and 400 of TEST600.GB: the 1001st
## pair generated (260/450) is row 7 of the best 20
DESCRIPTION="find-pair --report: the best 20 of the first 1001 pairs (W10131)"
"${PRIMERPICKER}" --config "${TRUN}/PD2_CRITERIA_W10131.CFG" find-pair --report "${REPORT}" \
    "${TEST600}" 260 400 > /dev/null && \
    [[ "$(body "${REPORT}"; printf x)" == "$(body "${TRUN}/TEST600_260_400_CRITERIA_W10131.PRN"; printf x)" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="find-pair: the list stops at 1000 pairs and warns"
"${PRIMERPICKER}" --config "${TRUN}/PD2_CRITERIA_W10131.CFG" find-pair "${TEST600}" 260 400 | \
    grep -qxF "More than 1000 possible pairs found!" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit 0
