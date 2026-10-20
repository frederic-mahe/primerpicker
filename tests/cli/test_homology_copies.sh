#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "primerpicker homology" on primers with copies
## elsewhere in the molecule: only the primer's own site is left out.
## Expected files come from PRIMER.exe (tests/legacy/reports/homology/;
## tests/legacy/README.md, "Copies of the primer elsewhere").
## Usage: test_homology_copies.sh /path/to/primerpicker /path/to/tests/legacy

## Print a header
SCRIPT_NAME="Test homology copies"
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

DESCRIPTION="homology --report: every recorded report of a primer with copies"
FAILED=""
while IFS=$'\t' read -r FILE MOLECULE POSITION KIND CONFIG ; do
    [[ "${FILE}" == \#* ]] && continue
    OPTIONS=()
    [[ "${CONFIG}" != "-" ]] && OPTIONS+=(--config "${LEGACY}/config/${CONFIG}")
    OPTIONS+=(homology --report "${REPORT}")
    [[ "${KIND}" == "false_priming" ]] && OPTIONS+=(--false-priming)
    "${PRIMERPICKER}" "${OPTIONS[@]}" "${LEGACY}/${MOLECULE}" "${POSITION}" > /dev/null && \
        same_body "${REPORT}" "${REPORTS}/${FILE}" || \
            FAILED="${FAILED} ${FILE}"
done < "${REPORTS}/homology_reports_copies.tsv"
[[ -z "${FAILED}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}:${FAILED}"

exit 0
