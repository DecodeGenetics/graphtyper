cmake_minimum_required(VERSION 3.2)

# Compile test.cpp as the main Catch file
add_library(catch-gt OBJECT test.cpp)

## Graph tests
set(graphtyper_graph_TEST_FILES
  graph/test_graph.cpp
  graph/test_constructor.cpp
  graph/test_genomic_region.cpp
  graph/test_haplotypes.cpp
)

add_executable(test_graphtyper_graph
  ${graphtyper_graph_TEST_FILES}
  $<TARGET_OBJECTS:catch-gt>
  $<TARGET_OBJECTS:graphtyper_objects>
)

target_link_libraries(test_graphtyper_graph ${graphtyper_all_libraries})

## Index tests
set(graphtyper_index_TEST_FILES
  index/test_index.cpp
)

add_executable(test_graphtyper_index
  ${graphtyper_index_TEST_FILES}
  $<TARGET_OBJECTS:catch-gt>
  $<TARGET_OBJECTS:graphtyper_objects>
)

target_link_libraries(test_graphtyper_index ${graphtyper_all_libraries})

## Typer tests
set(graphtyper_typer_TEST_FILES
  typer/test_path.cpp
  typer/test_genotype_path.cpp
  typer/test_vcf.cpp
  typer/test_vcf_io.cpp
)

add_executable(test_graphtyper_typer
  ${graphtyper_typer_TEST_FILES}
  $<TARGET_OBJECTS:catch-gt>
  $<TARGET_OBJECTS:graphtyper_objects>
)

target_link_libraries(test_graphtyper_typer ${graphtyper_all_libraries})

## Utilities tests
set(graphtyper_utilities_TEST_FILES
  utilities/test_kmer_help_functions.cpp
  utilities/test_utilities.cpp
)

add_executable(test_graphtyper_utilities
  ${graphtyper_utilities_TEST_FILES}
  $<TARGET_OBJECTS:catch-gt>
  $<TARGET_OBJECTS:graphtyper_objects>
)

target_link_libraries(test_graphtyper_utilities ${graphtyper_all_libraries})


add_test(NAME UtilityTests COMMAND test_graphtyper_utilities)
add_test(NAME GraphTests COMMAND test_graphtyper_graph)
add_test(NAME IndexTests COMMAND test_graphtyper_index)
add_test(NAME TyperTests COMMAND test_graphtyper_typer)
