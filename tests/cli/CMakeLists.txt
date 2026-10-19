# Black-box tests of the command-line front-end (bash scripts, in the
# style of swarm-tests and vsearch-tests). Not run when cross-compiling:
# the scripts need the host's bash.
if(CMAKE_CROSSCOMPILING)
  return()
endif()

find_program(BASH_EXECUTABLE NAMES bash REQUIRED)
add_test(NAME cli_analyze
  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_analyze.sh $<TARGET_FILE:primerpicker>)
add_test(NAME cli_pair
  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_pair.sh $<TARGET_FILE:primerpicker>)
add_test(NAME cli_find_pair
  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_find_pair.sh $<TARGET_FILE:primerpicker>
          ${PROJECT_SOURCE_DIR}/tests/legacy)
add_test(NAME cli_find_enzymes
  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_find_enzymes.sh $<TARGET_FILE:primerpicker>
          ${PROJECT_SOURCE_DIR}/tests/legacy)
add_test(NAME cli_multi
  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_multi.sh $<TARGET_FILE:primerpicker>)
add_test(NAME cli_options
  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_options.sh $<TARGET_FILE:primerpicker>
          ${PROJECT_SOURCE_DIR}/tests/legacy)
