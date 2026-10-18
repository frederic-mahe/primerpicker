#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "primerpicker pair". Expected values come from
## PRIMER.exe screens (tests/legacy/pairs/pair_screens.tsv).
## Usage: test_pair.sh /path/to/primerpicker

## Print a header
SCRIPT_NAME="Test pair"
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

DESCRIPTION="pair without two primers is an error"
"${PRIMERPICKER}" pair GTAAAACGACGGCCAGT > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="pair rejects letters other than A, C, G and T"
"${PRIMERPICKER}" pair GTAAAACGACGGCCAGT AACAGCTANGACCATG > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

## CCCTGTCATACCAATCTACCCC and TGTTATGCGCGTGGGG
PAIR=$("${PRIMERPICKER}" pair CCCTGTCATACCAATCTACCCC TGTTATGCGCGTGGGG)

check_line () {
    local DESCRIPTION="pair: ${1}"
    grep -qxE "${2}" <<< "${PAIR}" && \
        success "${DESCRIPTION}" || \
            failure "${DESCRIPTION}"
}

check_line "primer A" "Primer A: +5'-CCCTGTCATACCAATCTACCCC-3'"
check_line "primer B" "Primer B: +5'-TGTTATGCGCGTGGGG-3'"
check_line "length" "Length +22 +16 +6 base difference"
check_line "GC %" "% GC +54 +62! +8 % difference"
check_line "Tm" "Tm °C +73 +70 +3°C difference"
check_line "3' G/C runs" "Runs of G/C at 3' end +4! +4! +Runs found"
check_line "runs" "Runs of bases \(any\) +3 +2"
check_line "hairpins" "Hairpins +OK +OK"
check_line "3' dimers" "3' end matches +1 +1 +4!"
check_line "dimers" "Any adj homol bases +2 +4 +4"

## ATTCTTAAAGGCTACATTAA and CCCCAGTCAAGTAATAAAGA: equal lengths, no comment
DESCRIPTION="no comment when lengths are equal"
"${PRIMERPICKER}" pair ATTCTTAAAGGCTACATTAA CCCCAGTCAAGTAATAAAGA | grep -qxE "Length +20 +20" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit 0
