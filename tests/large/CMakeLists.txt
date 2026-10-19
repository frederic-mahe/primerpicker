# Opt-in tests of the largest molecules (option PRIMERPICKER_LARGE_TESTS)
find_program(BASH_EXECUTABLE NAMES bash REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(NAME large_molecule_size
  COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_molecule_size.sh $<TARGET_FILE:primerpicker>)
set_tests_properties(large_molecule_size PROPERTIES LABELS large TIMEOUT 1800)
