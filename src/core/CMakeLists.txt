add_library(primerpicker_core
  sequence/nucleotide.cpp)
target_include_directories(primerpicker_core PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_compile_options(primerpicker_core PRIVATE ${PRIMERPICKER_WARNING_FLAGS})
