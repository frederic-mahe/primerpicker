add_library(primerpicker_core
  thermo/legacy_tm.cpp)
target_include_directories(primerpicker_core PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_compile_options(primerpicker_core PRIVATE ${PRIMERPICKER_WARNING_FLAGS})
