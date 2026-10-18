# The command-line front-end: a library (also linked by the fuzz
# harnesses) and the program
add_library(primerpicker_cli
  cli.cpp)
target_link_libraries(primerpicker_cli PUBLIC primerpicker_core)
target_compile_options(primerpicker_cli PRIVATE ${PRIMERPICKER_WARNING_FLAGS})
target_compile_definitions(primerpicker_cli PRIVATE PRIMERPICKER_VERSION="${PROJECT_VERSION}")

add_executable(primerpicker
  main.cpp)
target_link_libraries(primerpicker PRIVATE primerpicker_cli)
target_compile_options(primerpicker PRIVATE ${PRIMERPICKER_WARNING_FLAGS})
