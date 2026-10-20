#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "--report FILE": the text reports of PRIMER.exe
## (Print to Disk), byte for byte but for the date. Expected files come
## from PRIMER.exe (tests/legacy/reports/; tests/legacy/README.md,
## "Reports").
## Usage: test_report.sh /path/to/primerpicker /path/to/tests/legacy


## Print a header
SCRIPT_NAME="Test reports"
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

## the report without its first line (the date)
body () {
    tail -n +2 "${1}"
}

## errors

DESCRIPTION="analyze --report needs a file name"
"${PRIMERPICKER}" analyze --report > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="analyze --report fails when the file cannot be written"
"${PRIMERPICKER}" analyze --report "${WORKDIR}/missing/REPORT" GTAAAACGACGGCCAGT > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="analyze --report given twice is an error"
"${PRIMERPICKER}" analyze --report "${REPORT}" --report "${REPORT}" GTAAAACGACGGCCAGT > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

## the One Primer report

DESCRIPTION="analyze --report writes INSPECT.PRN (but for the date)"
rm -f "${REPORT}"
"${PRIMERPICKER}" analyze --report "${REPORT}" GTAAAACGACGGCCAGT > /dev/null && \
    cmp -s <(body "${REPORT}") <(body "${REPORTS}/INSPECT.PRN") && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="analyze --report still prints the analysis"
"${PRIMERPICKER}" analyze --report "${REPORT}" GTAAAACGACGGCCAGT | \
    grep -q "^Sequence:      5'-GTAAAACGACGGCCAGT-3'$" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## the local date, or the UTC date where the standard library has no
## time zones (libc++)
DESCRIPTION="analyze --report: the first line holds today's date (MM/DD/YY)"
"${PRIMERPICKER}" analyze --report "${REPORT}" GTAAAACGACGGCCAGT > /dev/null
TITLE="Primer Design Report:  One Primer                     Report Date: "
FIRST_LINE=$(head -n 1 "${REPORT}")
[[ "${FIRST_LINE}" == "${TITLE}$(date +%m/%d/%y)"$'\r' || \
       "${FIRST_LINE}" == "${TITLE}$(date -u +%m/%d/%y)"$'\r' ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="analyze --report: options in any order"
rm -f "${REPORT}"
"${PRIMERPICKER}" analyze --extension "" --report "${REPORT}" GTAAAACGACGGCCAGT > /dev/null && \
    cmp -s <(body "${REPORT}") <(body "${REPORTS}/INSPECT.PRN") && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="analyze --report: every recorded Inspect report with pairs"
FAILED=""
while IFS=$'\t' read -r FILE PRIMER EXTENSION CONFIG DUPLEXES ; do
    [[ "${FILE}" == \#* || "${DUPLEXES}" == "no-pairs" ]] && continue
    OPTIONS=()
    [[ "${CONFIG}" != "-" ]] && OPTIONS+=(--config "${LEGACY}/config/${CONFIG}")
    OPTIONS+=(analyze --report "${REPORT}")
    [[ "${EXTENSION}" != "-" ]] && OPTIONS+=(--extension "${EXTENSION}")
    "${PRIMERPICKER}" "${OPTIONS[@]}" "${PRIMER}" > /dev/null && \
        cmp -s <(body "${REPORT}") <(body "${REPORTS}/${FILE}") || \
            FAILED="${FAILED} ${FILE}"
done < "${REPORTS}/inspect_reports.tsv"
[[ -z "${FAILED}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}:${FAILED}"

## a duplex without pairs: PRIMER.exe writes lines of tens of thousands
## of spaces; primerpicker draws the partner under the primer (decided
## on 2026-10-20)
DESCRIPTION="analyze --report: no long line for a primer without pairs"
"${PRIMERPICKER}" analyze --report "${REPORT}" AAAAAAAAAAAA > /dev/null && \
    (( $(wc -L < "${REPORT}") <= 80 )) && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## the Primer Pair report

DESCRIPTION="pair --report writes PAIR.PRN (but for the date)"
rm -f "${REPORT}"
"${PRIMERPICKER}" pair --report "${REPORT}" GTAAAACGACGGCCAGT AACAGCTATGACCATG > /dev/null && \
    cmp -s <(body "${REPORT}") <(body "${REPORTS}/PAIR.PRN") && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="pair --report fails when the file cannot be written"
"${PRIMERPICKER}" pair --report "${WORKDIR}/missing/REPORT" GTAAAACGACGGCCAGT AACAGCTATGACCATG \
    > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="pair --report: every recorded Primer Pair report with pairs"
FAILED=""
while IFS=$'\t' read -r FILE PRIMER_A PRIMER_B CONFIG DUPLEXES ; do
    [[ "${FILE}" == \#* || "${DUPLEXES}" == "no-pairs" ]] && continue
    OPTIONS=()
    [[ "${CONFIG}" != "-" ]] && OPTIONS+=(--config "${LEGACY}/config/${CONFIG}")
    "${PRIMERPICKER}" "${OPTIONS[@]}" pair --report "${REPORT}" "${PRIMER_A}" "${PRIMER_B}" > /dev/null && \
        cmp -s <(body "${REPORT}") <(body "${REPORTS}/${FILE}") || \
            FAILED="${FAILED} ${FILE}"
done < "${REPORTS}/pair_reports.tsv"
[[ -z "${FAILED}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}:${FAILED}"

## the Multiple Primers report

DESCRIPTION="multi --report writes MULTI.PRN (but for the date)"
rm -f "${REPORT}"
"${PRIMERPICKER}" multi --report "${REPORT}" CTCAAAGGCTCAGACATCAG GATCCATGCAATGCTAGCTA \
    ACGTTGCAGGTCAAGCTTGC TTGGCCAAGTCGATCGGATC > /dev/null && \
    cmp -s <(body "${REPORT}") <(body "${REPORTS}/MULTI.PRN") && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="multi --report: every recorded Multi Primer report"
FAILED=""
while IFS=$'\t' read -r FILE SELECTION FIELD_A FIELD_B FIELD_C FIELD_D CONFIG ; do
    [[ "${FILE}" == \#* ]] && continue
    OPTIONS=()
    [[ "${CONFIG}" != "-" ]] && OPTIONS+=(--config "${LEGACY}/config/${CONFIG}")
    "${PRIMERPICKER}" "${OPTIONS[@]}" multi --report "${REPORT}" \
        "${FIELD_A}" "${FIELD_B}" "${FIELD_C}" "${FIELD_D}" > /dev/null && \
        cmp -s <(body "${REPORT}") <(body "${REPORTS}/${FILE}") || \
            FAILED="${FAILED} ${FILE}"
done < "${REPORTS}/multi_reports.tsv"
[[ -z "${FAILED}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}:${FAILED}"

## the Find Pair report

DESCRIPTION="find-pair --report: every recorded Find Pair report"
FAILED=""
while IFS=$'\t' read -r FILE MOLECULE TARGET_A TARGET_B CONFIG ; do
    [[ "${FILE}" == \#* ]] && continue
    OPTIONS=()
    [[ "${CONFIG}" != "-" ]] && OPTIONS+=(--config "${LEGACY}/config/${CONFIG}")
    "${PRIMERPICKER}" "${OPTIONS[@]}" find-pair --report "${REPORT}" "${LEGACY}/${MOLECULE}" "${TARGET_A}" \
        "${TARGET_B}" > /dev/null && \
        cmp -s <(body "${REPORT}") <(body "${REPORTS}/${FILE}") || \
            FAILED="${FAILED} ${FILE}"
done < "${REPORTS}/find_pair_reports.tsv"
[[ -z "${FAILED}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}:${FAILED}"

## tests/legacy/search/find_pair_no_match.tsv: no matching primers
DESCRIPTION="find-pair --report writes no report when no pair is found"
rm -f "${REPORT}"
"${PRIMERPICKER}" find-pair --report "${REPORT}" "${LEGACY}/enzymes/TLONG.GB" 120 300 > /dev/null 2>&1
[[ ! -e "${REPORT}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit 0
