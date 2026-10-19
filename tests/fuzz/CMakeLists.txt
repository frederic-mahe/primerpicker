# libFuzzer harnesses for the readers of untrusted input (proposal
# §3.10). Built with the "fuzz" preset only; campaigns run with
# run_campaign.sh.
foreach(target IN ITEMS fuzz_genbank fuzz_fasta fuzz_molecule_file fuzz_pd2 fuzz_site_search fuzz_find_pair fuzz_multi_primer fuzz_pd2_cfg)
  add_executable(${target} ${target}.cpp)
  target_link_libraries(${target} PRIVATE primerpicker_core)
  target_compile_options(${target} PRIVATE ${PRIMERPICKER_WARNING_FLAGS})
  target_link_options(${target} PRIVATE -fsanitize=fuzzer)
endforeach()

add_executable(fuzz_cli fuzz_cli.cpp)
target_link_libraries(fuzz_cli PRIVATE primerpicker_cli)
target_compile_options(fuzz_cli PRIVATE ${PRIMERPICKER_WARNING_FLAGS})
target_compile_definitions(fuzz_cli PRIVATE PRIMERPICKER_LEGACY_DIR="${PROJECT_SOURCE_DIR}/tests/legacy")
target_link_options(fuzz_cli PRIVATE -fsanitize=fuzzer)
