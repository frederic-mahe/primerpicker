#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "find-pair --best N": the Find Pair reports of
## PRIMER.exe printed with a typed "Print best: [N]", byte for byte but for
## the date (tests/legacy/reports/find_pair_best/; tests/legacy/README.md,
## "Reports of more than 20 pairs").
## Usage: test_find_pair_best.sh /path/to/primerpicker /path/to/tests/legacy


## Print a header
SCRIPT_NAME="Test find-pair --best"
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
REPORTS="${LEGACY}/reports"
WORKDIR=$(mktemp -d)
trap 'rm -rf "${WORKDIR}"' EXIT
REPORT="${WORKDIR}/REPORT"

## the report without its dates: every page header has one ("Report
## Date: MM/DD/YY", on the first line and on the Primer Pair report's
## page 2)
body () {
    sed 's|Report Date: [0-9/]*|Report Date: DATE|' "${1}"
}

## whether two reports are the same but for their dates; bash
## only (no cmp: MSYS2 does not install diffutils by default). The
## final x keeps trailing line ends in the comparison
same_body () {
    [[ "$(body "${1}"; printf x)" == "$(body "${2}"; printf x)" ]]
}


BEST="${REPORTS}/find_pair_best"

DESCRIPTION="find-pair --best N: the recorded reports of 25 to 9999 pairs"
FAILED=""
for N in 25 30 50 9999 ; do
    "${PRIMERPICKER}" find-pair --best "${N}" --report "${REPORT}" "${LEGACY}/search/TFP.GB" 150 400 > /dev/null && \
        same_body "${REPORT}" "${BEST}/TFP_150_400_best${N}.PRN" || \
            FAILED="${FAILED} TFP_${N}"
done
for N in 100 1000 ; do
    "${PRIMERPICKER}" --config "${LEGACY}/config/PD2_CRITERIA.CFG" find-pair --best "${N}" --report "${REPORT}" \
        "${LEGACY}/pairs/TEST600.GB" 150 400 > /dev/null 2>&1 && \
        same_body "${REPORT}" "${BEST}/TEST600_150_400_CRITERIA_best${N}.PRN" || \
            FAILED="${FAILED} TEST600_${N}"
done
[[ -z "${FAILED}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}:${FAILED}"

DESCRIPTION="find-pair --best N: N from 1 to 9999"
"${PRIMERPICKER}" find-pair --best 0 --report "${REPORT}" "${LEGACY}/search/TFP.GB" 150 400 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}: 0 accepted"
"${PRIMERPICKER}" find-pair --best 10000 --report "${REPORT}" "${LEGACY}/search/TFP.GB" 150 400 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}: 10000 accepted"
"${PRIMERPICKER}" find-pair --best x --report "${REPORT}" "${LEGACY}/search/TFP.GB" 150 400 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}: x accepted"
success "${DESCRIPTION}"

DESCRIPTION="find-pair --best N needs --report"
"${PRIMERPICKER}" find-pair --best 25 "${LEGACY}/search/TFP.GB" 150 400 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}"
success "${DESCRIPTION}"

exit 0
