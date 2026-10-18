include(FetchContent)

# Use the system's Catch2 when there is one (find_package first), and
# download it otherwise. Presets that cannot link against the system
# library (libc++) set FETCHCONTENT_TRY_FIND_PACKAGE_MODE to NEVER.
FetchContent_Declare(Catch2
  URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.7.1.tar.gz
  URL_HASH SHA256=c991b247a1a0d7bb9c39aa35faf0fe9e19764213f28ffba3109388e62ee0269c
  SYSTEM
  FIND_PACKAGE_ARGS 3)
FetchContent_MakeAvailable(Catch2)

if(NOT COMMAND catch_discover_tests)
  # Catch.cmake is only on the module path when Catch2 was found with
  # find_package
  list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
endif()
include(Catch)

add_executable(unit_tests
  test_skeleton.cpp
  test_nucleotide.cpp
  test_breslauer_1986.cpp
  test_legacy_tm.cpp
  test_legacy_tm_primer3.cpp
  test_legacy_criteria.cpp
  test_genbank.cpp
  test_legacy_verdict.cpp
  test_legacy_hairpin.cpp
  test_legacy_pair.cpp
  test_pd2.cpp
  test_legacy_homology.cpp
  test_legacy_match.cpp
  test_legacy_find_pair.cpp
  test_legacy_find_pair_ranking.cpp
  test_rebase_enzymes.cpp
  test_iupac.cpp
  test_site_search.cpp
  test_molecule_topology.cpp
  test_molecule_formats.cpp)
target_link_libraries(unit_tests PRIVATE primerpicker_core Catch2::Catch2WithMain)
target_compile_options(unit_tests PRIVATE ${PRIMERPICKER_WARNING_FLAGS})
target_compile_definitions(unit_tests PRIVATE
  PRIMERPICKER_LEGACY_DIR="${PROJECT_SOURCE_DIR}/tests/legacy"
  PRIMERPICKER_TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/tests/data")

# PRE_TEST: list the tests when ctest runs, not after linking, so that
# a cross-compiled build succeeds even where the emulator (Wine)
# cannot run
catch_discover_tests(unit_tests DISCOVERY_MODE PRE_TEST)
