#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "primerpicker analyze". Expected values come from
## PRIMER.exe screens (tests/legacy/), unless stated otherwise.
## Usage: test_analyze.sh /path/to/primerpicker

## Print a header
SCRIPT_NAME="Test analyze"
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

DESCRIPTION="check if primerpicker is executable"
[[ -x "${PRIMERPICKER}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"


#*****************************************************************************#
#                                                                             #
#                             Command line                                    #
#                                                                             #
#*****************************************************************************#

DESCRIPTION="--version prints the version and succeeds"
"${PRIMERPICKER}" --version | grep -qE "^primerpicker [0-9]+\.[0-9]+\.[0-9]+$" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="--help prints the usage and succeeds"
"${PRIMERPICKER}" --help | grep -q "^Usage: primerpicker analyze \[--extension EXTENSION\] PRIMER" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="no argument is an error"
"${PRIMERPICKER}" > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="an unknown command is an error"
"${PRIMERPICKER}" foo > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="analyze without a primer is an error"
"${PRIMERPICKER}" analyze > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="analyze rejects letters other than A, C, G and T"
"${PRIMERPICKER}" analyze GTAAAACGNCGGCCAGT > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="analyze explains a rejected primer on stderr"
"${PRIMERPICKER}" analyze GTAAAACGNCGGCCAGT 2>&1 > /dev/null | grep -q "N" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="analyze rejects a primer of one base"
"${PRIMERPICKER}" analyze A > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="analyze rejects a primer longer than 50 bases (Direct Entry field)"
"${PRIMERPICKER}" analyze "$(printf 'A%.0s' {1..51})" > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"


#*****************************************************************************#
#                                                                             #
#                    Primer Analysis screen of PRIMER.exe                     #
#                                                                             #
#*****************************************************************************#

## M13 primer, Inspect screen
M13=$("${PRIMERPICKER}" analyze GTAAAACGACGGCCAGT)

check_line () {
    local DESCRIPTION="M13: ${1}"
    grep -qxE "${2}" <<< "${M13}" && \
        success "${DESCRIPTION}" || \
            failure "${DESCRIPTION}"
}

check_line "sequence" "Sequence: +5'-GTAAAACGACGGCCAGT-3'"
check_line "length" "Length: +17"
check_line "no 5' extension" "5' extension: +--none--"
check_line "composition" "Composition: +A:6 +C:4 +G:5 +T:2"
check_line "GC %" "% GC +min 50, max 60 +52 +YES"
check_line "Tm" "Tm °C +min 55, max 80 +64 +YES"
check_line "hairpins" "No hairpins +energy cutoff 0\.0 kcal +-- +YES"
check_line "3' dimers" "No 3' dimers +reject >= 3 matches at 3' end +2 +YES"
check_line "dimers" "No dimers +reject >= 7 adjacent homologous bases +4 +YES"
check_line "runs" "No runs +reject >= 3 base runs +4 +NO"
check_line "3' GC runs" "No 3' GC runs +reject >= 3 G or C at 3' end +0 +YES"
check_line "hairpin line" "Hairpin: +no hairpins found"
check_line "verdict" "Meets criteria: NO"

DESCRIPTION="lowercase primers are accepted (Direct Entry)"
"${PRIMERPICKER}" analyze gtaaaacgacggccagt | grep -qxE "Sequence: +5'-GTAAAACGACGGCCAGT-3'" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## Hairpin: GACGCTTTTGCGTC, -3.9 kcal (hairpin_design.tsv)
DESCRIPTION="the hairpin energy is shown"
"${PRIMERPICKER}" analyze GACGCTTTTGCGTC | grep -qxE "Hairpin: +-3\.9 kcal" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="a primer with a hairpin fails the hairpin criterion"
"${PRIMERPICKER}" analyze GACGCTTTTGCGTC | grep -qxE "No hairpins +energy cutoff 0\.0 kcal +-- +NO" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## A Tm that is not positive is shown as 0 (AAGATTTA, Inspect screen)
DESCRIPTION="a Tm that is not positive is shown as 0"
"${PRIMERPICKER}" analyze AAGATTTA | grep -qxE "Tm °C +min 55, max 80 +0 +NO" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## A primer that meets every criterion: PRIMER.exe shows "MEETS CRITERIA
## SELECTED: YES" (recorded with inspect_primer.sh; its values are in
## tm/direct_entry_tm.tsv and inspect/inspect.tsv)
DESCRIPTION="a primer can meet all criteria"
"${PRIMERPICKER}" analyze ACACACACACACACACACAC | grep -qx "Meets criteria: YES" && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"


#*****************************************************************************#
#                                                                             #
#                               5' extension                                  #
#                                                                             #
#*****************************************************************************#

## CCCC + GTAAAACGACGGCCAGT, Inspect screen (inspect_extension.tsv)
EXTENDED=$("${PRIMERPICKER}" analyze --extension CCCC GTAAAACGACGGCCAGT)

check_extended () {
    local DESCRIPTION="extension: ${1}"
    grep -qxE "${2}" <<< "${EXTENDED}" && \
        success "${DESCRIPTION}" || \
            failure "${DESCRIPTION}"
}

check_extended "sequence" "Sequence: +5'-CCCCGTAAAACGACGGCCAGT-3'"
check_extended "length includes the extension" "Length: +21"
check_extended "extension" "5' extension: +CCCC"
check_extended "GC % includes the extension" "% GC +min 50, max 60 +61 +NO"
check_extended "Tm without the extension" "Tm °C +min 55, max 80 +64 +YES"
check_extended "hairpin with the extension" "Hairpin: +-2\.7 kcal"
check_extended "runs" "No runs +reject >= 3 base runs +4 +NO"

DESCRIPTION="the extension is checked like a primer"
"${PRIMERPICKER}" analyze --extension CCNC GTAAAACGACGGCCAGT > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

exit 0
