cmake_minimum_required(VERSION 2.8.8)

set(graphtyper_index_TEST_FILES
  ../test.cpp
  test_index.cpp
)

add_executable(test_graphtyper_index ${graphtyper_index_TEST_FILES} $<TARGET_OBJECTS:graphtyper_objects>)
target_link_libraries(test_graphtyper_index ${graphtyper_all_libraries})
