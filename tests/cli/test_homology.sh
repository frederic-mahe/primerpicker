#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "primerpicker homology": the Uniqueness and False
## Priming screens of PRIMER.exe, and their reports. Expected values come
## from PRIMER.exe (tests/legacy/homology/homology.tsv, tests/legacy/reports/;
## tests/legacy/README.md, "Uniqueness and false priming").
## Usage: test_homology.sh /path/to/primerpicker /path/to/tests/legacy

## Print a header
SCRIPT_NAME="Test homology"
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

DESCRIPTION="homology needs a molecule and a position"
"${PRIMERPICKER}" homology "${TREP2}" > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="homology rejects position 0"
"${PRIMERPICKER}" homology "${TREP2}" 0 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

## TREP2 has 600 bases: a 20-base primer starts at 581 at most
DESCRIPTION="homology accepts the last position of a primer"
"${PRIMERPICKER}" homology "${TREP2}" 581 > /dev/null 2>&1 && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="homology rejects a primer beyond the molecule's end"
"${PRIMERPICKER}" homology "${TREP2}" 582 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="homology rejects a missing molecule file"
"${PRIMERPICKER}" homology "${WORKDIR}/missing.gb" 51 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

## the screens (tests/legacy/homology/homology.tsv: TREP2.GB, 51)

DESCRIPTION="homology lists the regions of homology in the screen's order"
[[ $("${PRIMERPICKER}" homology "${TREP2}" 51 | grep -E '^[ACGT]+ ' | tr -s ' ' | tr '\n' ';') == \
   "AAAGGCTCAGA 330 (C) 32;CAAAGGCTCA 270 27;TCAGACATCAG 230 (C) 18;CAGACATCAG 150 10;" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="homology --false-priming lists the false priming sites"
[[ $("${PRIMERPICKER}" homology --false-priming "${TREP2}" 51 | grep -E '^[ACGT]+ ' | tr -s ' ' | tr '\n' ';') == \
   "TCAGACATCAG 230 (C) 18;CAGACATCAG 150 10;" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## TEST600.GB, 51: "Primer is unique" (README, "Uniqueness and false
## priming")
DESCRIPTION="homology: a unique primer"
"${PRIMERPICKER}" homology "${LEGACY}/pairs/TEST600.GB" 51 | grep -qx "Primer is unique" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## the reports

DESCRIPTION="homology --report writes UNIQUE.PRN (but for the date)"
rm -f "${REPORT}"
"${PRIMERPICKER}" homology --report "${REPORT}" "${TREP2}" 51 > /dev/null && \
    same_body "${REPORT}" "${REPORTS}/UNIQUE.PRN" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="homology --false-priming --report writes FALSEPRI.PRN (but for the date)"
rm -f "${REPORT}"
"${PRIMERPICKER}" homology --false-priming --report "${REPORT}" "${TREP2}" 51 > /dev/null && \
    same_body "${REPORT}" "${REPORTS}/FALSEPRI.PRN" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="homology --report: every recorded report"
FAILED=""
while IFS=$'\t' read -r FILE MOLECULE POSITION KIND CONFIG NOTE ; do
    [[ "${FILE}" == \#* || "${NOTE}" != "-" ]] && continue
    OPTIONS=()
    [[ "${CONFIG}" != "-" ]] && OPTIONS+=(--config "${LEGACY}/config/${CONFIG}")
    OPTIONS+=(homology --report "${REPORT}")
    [[ "${KIND}" == "false_priming" ]] && OPTIONS+=(--false-priming)
    "${PRIMERPICKER}" "${OPTIONS[@]}" "${LEGACY}/${MOLECULE}" "${POSITION}" > /dev/null && \
        same_body "${REPORT}" "${REPORTS}/${FILE}" || \
            FAILED="${FAILED} ${FILE}"
done < "${REPORTS}/homology_reports.tsv"
[[ -z "${FAILED}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}:${FAILED}"

## PRIMER.exe shows a box instead ("No false priming found") and prints
## nothing
DESCRIPTION="homology --report writes no report when there is no site"
rm -f "${REPORT}"
"${PRIMERPICKER}" homology --report "${REPORT}" "${LEGACY}/pairs/TEST600.GB" 51 > /dev/null 2>&1
[[ ! -e "${REPORT}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit 0
