cmake_minimum_required(VERSION 2.8.8)

# Compile test.cpp as the main Catch file
add_library(catch OBJECT test.cpp)

add_subdirectory(graph)
add_subdirectory(index)
add_subdirectory(typer)
add_subdirectory(utilities)

add_test(NAME UtilityTests COMMAND test_graphtyper_utilities)
add_test(NAME GraphTests COMMAND test_graphtyper_graph)
add_test(NAME IndexTests COMMAND test_graphtyper_index)
add_test(NAME TyperTests COMMAND test_graphtyper_typer)
