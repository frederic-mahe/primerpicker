#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "primerpicker match": the Find Matching Primer list
## of PRIMER.exe (Match Active Primer), and its report. Expected values
## come from PRIMER.exe (tests/legacy/search/match_active.tsv,
## tests/legacy/search/match_3prime_runs.tsv, tests/legacy/reports/match/;
## tests/legacy/README.md, "Match Active Primer").
## Usage: test_match.sh /path/to/primerpicker /path/to/tests/legacy

## Print a header
SCRIPT_NAME="Test match"
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
TREP2="${LEGACY}/homology/TREP2.GB"
WORKDIR=$(mktemp -d)
trap 'rm -rf "${WORKDIR}"' EXIT
REPORT="${WORKDIR}/REPORT"

## the report without its first line (the date)
body () {
    tail -n +2 "${1}"
}

## whether two reports are the same but for their first line (bash
## only; the final x keeps trailing line ends in the comparison)
same_body () {
    [[ "$(body "${1}"; printf x)" == "$(body "${2}"; printf x)" ]]
}

## errors

DESCRIPTION="match needs a molecule, a position and a target"
"${PRIMERPICKER}" match "${LEGACY}/pairs/TEST600.GB" 51 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="match rejects a primer A beyond the molecule's end"
"${PRIMERPICKER}" match "${LEGACY}/pairs/TEST600.GB" 582 400 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="match rejects target 0"
"${PRIMERPICKER}" match "${LEGACY}/pairs/TEST600.GB" 51 0 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

## the lists (search/match_active.tsv: TEST600.GB, primer A at 51,
## target 400: 11 primers B)

DESCRIPTION="match lists the primers B of TEST600.GB, primer A at 51, target 400"
[[ $("${PRIMERPICKER}" match "${LEGACY}/pairs/TEST600.GB" 51 400 | awk '$1 ~ /^[0-9]+$/ && NF >= 8 {print $1}' | tr '\n' ',') == \
   "$(grep -v '^#' "${LEGACY}/search/match_active.tsv" | cut -f9 | tr '\n' ',')" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="match: the recorded lists of TRUN.GB (3'-end runs)"
FAILED=""
while IFS=$'\t' read -r MOLECULE POSITION TARGET LISTED ; do
    [[ "${MOLECULE}" == \#* ]] && continue
    FOUND=$("${PRIMERPICKER}" match "${LEGACY}/${MOLECULE}" "${POSITION}" "${TARGET}" | \
                awk '$1 ~ /^[0-9]+$/ && NF >= 8 {print $1}' | paste -sd, -)
    [[ "${FOUND}" == "${LISTED}" ]] || FAILED="${FAILED} ${TARGET}"
done < "${LEGACY}/search/match_3prime_runs.tsv"
[[ -z "${FAILED}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}:${FAILED}"

## the reports

DESCRIPTION="match --report: every recorded Match Active Primer report"
FAILED=""
while IFS=$'\t' read -r FILE MOLECULE POSITION TARGET CONFIG ; do
    [[ "${FILE}" == \#* ]] && continue
    OPTIONS=()
    [[ "${CONFIG}" != "-" ]] && OPTIONS+=(--config "${LEGACY}/config/${CONFIG}")
    "${PRIMERPICKER}" "${OPTIONS[@]}" match --report "${REPORT}" "${LEGACY}/${MOLECULE}" "${POSITION}" "${TARGET}" \
        > /dev/null && \
        same_body "${REPORT}" "${REPORTS}/${FILE}" || \
            FAILED="${FAILED} ${FILE}"
done < "${REPORTS}/match_reports.tsv"
[[ -z "${FAILED}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}:${FAILED}"

exit 0
