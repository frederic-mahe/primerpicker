add_executable(primerpicker
  main.cpp)
target_link_libraries(primerpicker PRIVATE primerpicker_core)
target_compile_options(primerpicker PRIVATE ${PRIMERPICKER_WARNING_FLAGS})
target_compile_definitions(primerpicker PRIVATE PRIMERPICKER_VERSION="${PROJECT_VERSION}")
