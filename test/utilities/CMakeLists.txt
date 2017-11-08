cmake_minimum_required(VERSION 2.8.8)

set(graphtyper_utilities_TEST_FILES
  test_adapter_removal.cpp
  test_kmer_help_functions.cpp
  test_utilities.cpp
)

add_executable(test_graphtyper_utilities
  ${graphtyper_utilities_TEST_FILES}
  $<TARGET_OBJECTS:catch>
  $<TARGET_OBJECTS:graphtyper_objects>
)

target_link_libraries(test_graphtyper_utilities ${graphtyper_all_libraries})
