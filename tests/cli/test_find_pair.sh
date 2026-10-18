#!/bin/bash -
# shellcheck disable=SC2015

## Black-box tests of "primerpicker find-pair". Expected values come
## from PRIMER.exe screens (tests/legacy/search/find_pair.tsv,
## find_pair_weights.tsv and find_pair_false_priming.tsv).
## Usage: test_find_pair.sh /path/to/primerpicker /path/to/tests/legacy

## Print a header
SCRIPT_NAME="Test find-pair"
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
TFP="${LEGACY}/search/TFP.GB"

DESCRIPTION="find-pair without targets is an error"
"${PRIMERPICKER}" find-pair "${TEST600}" 150 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="find-pair rejects a missing molecule file"
"${PRIMERPICKER}" find-pair "${LEGACY}/no_such_file.gb" 150 400 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="find-pair rejects a target that is not a number"
"${PRIMERPICKER}" find-pair "${TEST600}" 15O 400 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="find-pair rejects weights above 5"
"${PRIMERPICKER}" find-pair --weights 11161 "${TEST600}" 150 400 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

DESCRIPTION="find-pair rejects fewer than five weights"
"${PRIMERPICKER}" find-pair --weights 1111 "${TEST600}" 150 400 > /dev/null 2>&1 && \
    failure "${DESCRIPTION}" || \
        success "${DESCRIPTION}"

## TEST600.GB, targets 150 and 400, built-in weights
LIST=$("${PRIMERPICKER}" find-pair "${TEST600}" 150 400)

check_line () {
    local DESCRIPTION="find-pair: ${1}"
    grep -qxE "${2}" <<< "${LIST}" && \
        success "${DESCRIPTION}" || \
            failure "${DESCRIPTION}"
}

check_line "possible primers A" "Primer A: +3 possible primers"
check_line "possible primers B" "Primer B: +7 possible primers"
check_line "possible pairs" "Pairs: +10 possible"
check_line "first pair" " +73 +425 +50 / 50 +69 / 67 +2 +5 +- / -"
check_line "pair with a false-priming site" " +71 +453 +55 / 50 +75 / 69 +2 +4 +21 / -"

## the order of the list: positions A and B of the pairs, by rank
order () {
    "${PRIMERPICKER}" find-pair "$@" | \
        awk '/-A-/ {list = 1; next} list && NF > 0 {printf "%s/%s ", $1, $2}'
}

DESCRIPTION="find-pair: order with the built-in weights"
[[ $(order "${TEST600}" 150 400) == "73/425 73/424 73/413 73/414 72/425 72/424 72/413 72/414 71/453 71/457 " ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

DESCRIPTION="find-pair: order with weights 0 5 0 0 1"
[[ $(order --weights 05001 "${TEST600}" 150 400) == "73/425 73/424 72/425 72/424 73/413 73/414 72/413 72/414 71/453 71/457 " ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

## TFP.GB, weights 0 5 0 5 2: the ten last pairs of the list
DESCRIPTION="find-pair: order with false-priming sites (weights 0 5 0 5 2)"
TAIL=$(order --weights 05052 "${TFP}" 150 400 | awk '{for (i = NF - 9; i <= NF; i++) printf "%s ", $i}')
EXPECTED=$(awk -F'\t' '$4 == "0 5 0 5 2" && $5 > 20 {printf "%s/%s ", $6, $7}' "${LEGACY}/search/find_pair_false_priming.tsv")
[[ -n "${EXPECTED}" && "${TAIL}" == "${EXPECTED}" ]] && \
    success "${DESCRIPTION}" || \
        failure "${DESCRIPTION}"

exit 0
