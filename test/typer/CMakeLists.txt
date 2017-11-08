cmake_minimum_required(VERSION 2.8.8)

set(graphtyper_typer_TEST_FILES
  test_path.cpp
  test_genotype_path.cpp
  test_graph_swapping.cpp
  test_vcf.cpp
  test_vcf_io.cpp
)

add_executable(test_graphtyper_typer ${graphtyper_typer_TEST_FILES} $<TARGET_OBJECTS:catch> $<TARGET_OBJECTS:graphtyper_objects>)
target_link_libraries(test_graphtyper_typer ${graphtyper_all_libraries})
