#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "primerpicker explain": the built-in definitions
## (data/glossary/glossary.md, reviewed by the lab on 2026-10-21).
## Usage: test_explain.sh /path/to/primerpicker /path/to/data/glossary/glossary.md

## Print a header
SCRIPT_NAME="Test explain"
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
GLOSSARY="${2}"

DESCRIPTION="explain lists the 25 terms, one per line"
[[ $("${PRIMERPICKER}" explain | wc -l) -eq 25 ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="explain lists each term's name and title"
"${PRIMERPICKER}" explain | grep -qxE "gc +%GC \(GC content\)" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="explain TERM prints the entry's title and text"
OUTPUT=$("${PRIMERPICKER}" explain gc)
[[ "$(head -n 1 <<< "${OUTPUT}")" == "%GC (GC content)" ]] && \
    grep -qF "    %GC = floor(100 × (number of G + number of C) / length)" <<< "${OUTPUT}" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="explain accepts aliases, in any case"
"${PRIMERPICKER}" explain Ta | head -n 1 | grep -qxF "Annealing temperature (Ta)" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="explain TERM prints the text of glossary.md as is"
EXPECTED=$(awk '/^## Hairpin$/ {found = 1; next} /^## / {found = 0} found' "${GLOSSARY}" | sed -e '/./,$!d')
ACTUAL=$("${PRIMERPICKER}" explain hairpin | tail -n +3)
[[ "${ACTUAL}" == "$(printf '%s' "${EXPECTED}")" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="explain rejects an unknown term"
"${PRIMERPICKER}" explain no-such-term > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="explain takes at most one term"
"${PRIMERPICKER}" explain gc tm > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

exit 0
