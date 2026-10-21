add_library(primerpicker_core
  colour/colour.cpp
  report/legacy_report.cpp
  thermo/legacy_tm.cpp)
target_include_directories(primerpicker_core PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_compile_options(primerpicker_core PRIVATE ${PRIMERPICKER_WARNING_FLAGS})
