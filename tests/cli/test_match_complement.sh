#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "match --complement" (primer A on the complement
## strand). Expected values come from PRIMER.exe
## (tests/legacy/complement/*.PRN, recorded 2026-10-22;
## tests/legacy/PROBES_20261022_complement.md).
## Usage: test_match_complement.sh /path/to/primerpicker /path/to/tests/legacy

## Print a header
SCRIPT_NAME="Test match --complement"
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
MATCH_ALL="${LEGACY}/match_variations/PD2_MATCH_ALL.CFG"
WORKDIR=$(mktemp -d)
trap 'rm -rf "${WORKDIR}"' EXIT
REPORT="${WORKDIR}/REPORT"

## the report without its date (bash only: no cmp in MSYS2)
body () {
    sed 's|Report Date: [0-9/]*|Report Date: DATE|' "${1}"
}

## primer A at 400(C) and 460(C), primers B on the top strand (the
## form's default for such a primer A), PD2_MATCH_ALL.CFG
for RUN in "400 100 match_c400_all" "460 300 match_c460_all" ; do
    read -r POSITION TARGET NAME <<< "${RUN}"
    DESCRIPTION="match --complement --report: PRIMER.exe's report ${NAME}, but for the date"
    "${PRIMERPICKER}" --config "${MATCH_ALL}" match --complement --report "${REPORT}" \
                      "${TEST600}" "${POSITION}" "${TARGET}" > /dev/null && \
        [[ "$(body "${REPORT}"; printf x)" == "$(body "${LEGACY}/complement/${NAME}.PRN"; printf x)" ]] && \
        success "${DESCRIPTION}" || \
            failure "${DESCRIPTION}"
done

## built-in settings: primer A's row has no '!' (GC 45 at 400(C))
DESCRIPTION="match --complement --report: PRIMER.exe's report match_c400, but for the date"
"${PRIMERPICKER}" match --complement --report "${REPORT}" "${TEST600}" 400 100 > /dev/null && \
    [[ "$(body "${REPORT}"; printf x)" == "$(body "${LEGACY}/complement/match_c400.PRN"; printf x)" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="match --complement: primers B on the top strand by default"
"${PRIMERPICKER}" match --complement "${TEST600}" 400 100 | \
    grep -q "primers B on the top strand" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="match --complement: the screen names primer A's strand"
"${PRIMERPICKER}" match --complement "${TEST600}" 400 100 | \
    grep -q "position 400 (complement strand)" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## match_c400.PRN: one primer B, 73
DESCRIPTION="match --complement: the recorded primer B"
"${PRIMERPICKER}" match --complement "${TEST600}" 400 100 | \
    grep -qE "^ +73 +27 +50 +69 +2 +4 +1 +3 +2 +0$" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="match: --complement given twice is an error"
"${PRIMERPICKER}" match --complement --complement "${TEST600}" 400 100 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="match --complement: a primer A that does not fit is an error"
"${PRIMERPICKER}" match --complement "${TEST600}" 19 100 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

exit 0
