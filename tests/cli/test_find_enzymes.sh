#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "primerpicker find-enzymes". Expected values come
## from PRIMER.exe screens (tests/legacy/enzymes/find_enzymes.tsv).
## Usage: test_find_enzymes.sh /path/to/primerpicker /path/to/tests/legacy

## Print a header
SCRIPT_NAME="Test find-enzymes"
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
TENZ="${LEGACY}/enzymes/TENZ.GB"
RECORDED="${LEGACY}/enzymes/find_enzymes.tsv"

DESCRIPTION="find-enzymes without an enzyme is an error"
"${PRIMERPICKER}" find-enzymes "${TENZ}" > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="find-enzymes rejects an unknown enzyme"
"${PRIMERPICKER}" find-enzymes "${TENZ}" NoSuchI > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="find-enzymes rejects a missing molecule file"
"${PRIMERPICKER}" find-enzymes "${LEGACY}/no_such_file.gb" AluI > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="find-enzymes: header of the true sites"
"${PRIMERPICKER}" find-enzymes "${TENZ}" AlwI | \
    grep -qxE "Search for: AlwI +GGATC +Number of Sites = 2" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="find-enzymes: header of the almost sites"
"${PRIMERPICKER}" find-enzymes --almost "${TENZ}" AlwI | \
    grep -qxE "Search for: AlwI +GGATC +\(GATCC\) +Number of Sites = 18" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## the sites as find_enzymes.tsv writes them: position, strand, match
sites () {
    "${PRIMERPICKER}" find-enzymes "$@" | \
        awk '/^Search for:/ {list = 1; next}
             list && /\(/ {gsub(/[()]/, "", $2); print $1, "complement", $2; next}
             list && /C$/ {sub(/C$/, "", $1); print $1, "complement", "-"; next}
             list && NF == 2 {print $1, "top", $2; next}
             list {print $1, "top", "-"}'
}

recorded () {
    awk -F'\t' -v kind="$1" -v enzyme="$2" \
        '$1 == "TENZ.GB" && $2 == kind && $3 == enzyme {print $6, $7, $8}' "${RECORDED}"
}

for ENZYME in AluI AlwI FokI ApaI ; do
    DESCRIPTION="find-enzymes: true sites of ${ENZYME}"
    EXPECTED=$(recorded true "${ENZYME}")
    [[ -n "${EXPECTED}" && $(sites "${TENZ}" "${ENZYME}") == "${EXPECTED}" ]] && \
        success "${DESCRIPTION}" || \
            failure "${DESCRIPTION}"
done

for ENZYME in AlwI ApaI ; do
    DESCRIPTION="find-enzymes: almost sites of ${ENZYME}"
    EXPECTED=$(recorded almost "${ENZYME}")
    [[ -n "${EXPECTED}" && $(sites --almost "${TENZ}" "${ENZYME}") == "${EXPECTED}" ]] && \
        success "${DESCRIPTION}" || \
            failure "${DESCRIPTION}"
done

## circular molecule: a site across the origin, and sites on both
## strands at the same position
TENZC="${LEGACY}/enzymes/TENZC.GB"
recorded_on () {
    awk -F'\t' -v molecule="$1" -v kind="$2" -v enzyme="$3" \
        '$1 == molecule && $2 == kind && $3 == enzyme {print $6, $7, $8}' "${RECORDED}"
}

DESCRIPTION="find-enzymes: true sites across the origin of a circular molecule"
EXPECTED=$(recorded_on TENZC.GB true AluI)
[[ -n "${EXPECTED}" && $(sites "${TENZC}" AluI) == "${EXPECTED}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="find-enzymes: almost sites at the same position on both strands"
EXPECTED=$(recorded_on TENZC.GB almost PleI)
[[ -n "${EXPECTED}" && $(sites --almost "${TENZC}" PleI) == "${EXPECTED}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit 0
