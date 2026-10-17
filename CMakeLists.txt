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

if(BUILD_TESTING)
  add_subdirectory(tests/unit)
endif()
