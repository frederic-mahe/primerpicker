cmake_minimum_required(VERSION 3.28)

project(primerpicker
  VERSION 0.1.0
  DESCRIPTION "Interactive PCR and sequencing primer design"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# no C++20 modules: scanning would need an unversioned clang-scan-deps
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

# Warning flags are chosen per compiler in CMakePresets.json, and are
# applied to the project's own targets only (not to dependencies)
set(PRIMERPICKER_WARNING_FLAGS "" CACHE STRING
  "Warning flags for primerpicker targets (semicolon-separated list)")

include(CTest)

add_subdirectory(src/core)
add_subdirectory(src/cli)

if(BUILD_TESTING)
  add_subdirectory(tests/unit)
  add_subdirectory(tests/cli)
endif()

# Tests of the largest molecules (more than 2^32 bases: about 13 GB of
# memory and 4.4 GB of disk); run with ctest -L large
option(PRIMERPICKER_LARGE_TESTS "Add the tests of the largest molecules" OFF)
if(BUILD_TESTING AND PRIMERPICKER_LARGE_TESTS)
  add_subdirectory(tests/large)
endif()

option(PRIMERPICKER_FUZZ "Build the libFuzzer harnesses (clang only)" OFF)
if(PRIMERPICKER_FUZZ)
  add_subdirectory(tests/fuzz)
endif()
