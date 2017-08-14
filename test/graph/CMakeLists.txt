cmake_minimum_required(VERSION 2.8.8)

set(graphtyper_graph_TEST_FILES
  test_graph.cpp
  test_constructor.cpp
  test_genomic_region.cpp
  test_haplotypes.cpp
)

add_executable(test_graphtyper_graph ${graphtyper_graph_TEST_FILES} $<TARGET_OBJECTS:catch> $<TARGET_OBJECTS:graphtyper_objects>)
target_link_libraries(test_graphtyper_graph ${graphtyper_all_libraries})
