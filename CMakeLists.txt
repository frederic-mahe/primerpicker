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

# The graphical front-end (Qt 6 Widgets, 6.8 or later for the colour
# scheme setting): AUTO builds it when Qt is found, ON requires it
set(PRIMERPICKER_GUI "AUTO" CACHE STRING "Build the graphical front-end: AUTO, ON or OFF")
set_property(CACHE PRIMERPICKER_GUI PROPERTY STRINGS AUTO ON OFF)
set(PRIMERPICKER_QT_MINIMUM_VERSION 6.8)
set(primerpicker_qt_components Widgets)
if(BUILD_TESTING)
  list(APPEND primerpicker_qt_components Test)
endif()
if(PRIMERPICKER_GUI STREQUAL "ON")
  find_package(Qt6 ${PRIMERPICKER_QT_MINIMUM_VERSION} REQUIRED COMPONENTS ${primerpicker_qt_components})
elseif(PRIMERPICKER_GUI STREQUAL "AUTO")
  find_package(Qt6 ${PRIMERPICKER_QT_MINIMUM_VERSION} QUIET COMPONENTS ${primerpicker_qt_components})
elseif(NOT PRIMERPICKER_GUI STREQUAL "OFF")
  message(FATAL_ERROR "PRIMERPICKER_GUI must be AUTO, ON or OFF, not '${PRIMERPICKER_GUI}'")
endif()
if(Qt6_FOUND)
  message(STATUS "Graphical front-end: built with Qt ${Qt6_VERSION}")
else()
  message(STATUS "Graphical front-end: not built (PRIMERPICKER_GUI=${PRIMERPICKER_GUI})")
endif()

add_subdirectory(src/core)
add_subdirectory(src/cli)
if(Qt6_FOUND)
  add_subdirectory(src/gui)
endif()

if(BUILD_TESTING)
  add_subdirectory(tests/unit)
  add_subdirectory(tests/cli)
  if(Qt6_FOUND)
    add_subdirectory(tests/gui)
  endif()
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
