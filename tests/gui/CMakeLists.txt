# Tests of the graphical front-end: Catch2 test cases, with the input
# simulation of Qt Test, run without a display (offscreen
# platform, chosen by main.cpp)

# Catch2 is set up by tests/unit. A system Catch2 is found again here
# (imported targets are local to a directory), and puts its CMake
# module (Catch.cmake) on the module path; a downloaded one is global.
if(NOT TARGET Catch2::Catch2)
  find_package(Catch2 3 REQUIRED)
endif()
if(NOT COMMAND catch_discover_tests)
  FetchContent_GetProperties(Catch2)
  list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
endif()
include(Catch)

add_executable(gui_tests
  main.cpp
  test_main_window.cpp
  test_direct_entry.cpp
  test_from_molecule.cpp
  test_settings.cpp)
target_link_libraries(gui_tests PRIVATE primerpicker_gui Qt6::Test Catch2::Catch2)
target_compile_options(gui_tests PRIVATE ${PRIMERPICKER_WARNING_FLAGS})
target_compile_definitions(gui_tests PRIVATE PRIMERPICKER_LEGACY_DIR="${PROJECT_SOURCE_DIR}/tests/legacy")

catch_discover_tests(gui_tests DISCOVERY_MODE PRE_TEST TEST_PREFIX "gui: ")
