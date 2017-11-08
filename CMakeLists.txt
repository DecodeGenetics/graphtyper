cmake_minimum_required (VERSION 2.8.11)
project (graphtyper)

include(ExternalProject)

# The version number.
set (graphtyper_VERSION_MAJOR 1)
set (graphtyper_VERSION_MINOR 2)

# SeqAn options
set(SEQAN_BUILD_SYSTEM "SEQAN_RELEASE_LIBRARY" CACHE TYPE STRING)

if ($ENV{GT_BUILD_TYPE})
  set(GT_BUILD_TYPE $ENV{GT_BUILD_TYPE})
endif()

message(STATUS "Build type: ${GT_BUILD_TYPE}")
set(GT_COMMON_FLAGS "-Wall -Wextra -Wfatal-errors -pedantic -Wno-variadic-macros -std=c++11 -DSEQAN_HAS_ZLIB=1 -DSEQAN_USE_HTSLIB=1 -DSEQAN_ENABLE_TESTING=0 -march=core2 -mtune=generic")

if(GT_BUILD_TYPE STREQUAL "DEBUG" OR CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Building in debug mode.")
  set (CMAKE_CXX_FLAGS "${GT_COMMON_FLAGS} -g -O0 -DDEBUG -DSEQAN_ENABLE_DEBUG=1")
elseif(GT_BUILD_TYPE STREQUAL "RELDEB" OR GT_BUILD_TYPE STREQUAL "Reldeb")
  message(STATUS "Building in release (with assertions and debug info) mode.")
  set (CMAKE_CXX_FLAGS "${GT_COMMON_FLAGS} -g -O3 -DDEBUG -DSEQAN_ENABLE_DEBUG=1")
elseif(GT_BUILD_TYPE STREQUAL "PROFILING" OR GT_BUILD_TYPE STREQUAL "Profiling")
  message(STATUS "Building in profiling mode, which is release mode with profiling mode enabled.")
  set (CMAKE_CXX_FLAGS "${GT_COMMON_FLAGS} -pg -O3 -DNDEBUG -DSEQAN_ENABLE_DEBUG=0")
else()
  message(STATUS "Building in release mode.")
  set (CMAKE_CXX_FLAGS "${GT_COMMON_FLAGS} -O3 -DNDEBUG -DSEQAN_ENABLE_DEBUG=0")
endif()

message(STATUS "CXX flags are: ${CMAKE_CXX_FLAGS}")

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/utilities/cmake_modules/")

# configure a header file to pass some of the CMake settings to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/graphtyper/constants.hpp.in"
  "${PROJECT_BINARY_DIR}/graphtyper/constants.hpp"
)

include_directories("${PROJECT_BINARY_DIR}")

# Find external libraries
# zlib
message (STATUS "Checking for zlib")
find_package(ZLIB REQUIRED)
include_directories(SYSTEM ${ZLIB_INCLUDE_DIRS})

# bzip2
message (STATUS "Checking for bzip2")
find_package (BZip2 REQUIRED)
include_directories(SYSTEM ${BZIP_INCLUDE_DIRS})

# boost
set(Boost_USE_STATIC_LIBS ON)
message (STATUS "Checking for Boost")
find_package(Boost 1.57.0 REQUIRED COMPONENTS system iostreams log thread serialization filesystem REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# snappy
ExternalProject_Add(
  project_snappy
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/snappy
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/snappy
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/snappy/autogen.sh COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/snappy/configure
  BUILD_COMMAND $(MAKE)
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(project_snappy install_dir)
add_library(snappy STATIC IMPORTED)
include_directories(SYSTEM ${install_dir})

set(snappy_location ${install_dir}/.libs/libsnappy.a)
message(STATUS "snappy target location is ${snappy_location}")
set_property(TARGET snappy PROPERTY IMPORTED_LOCATION ${snappy_location})
add_dependencies(snappy project_snappy)

## rocksdb

ExternalProject_Add(
  project_rocksdb
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/rocksdb
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/rocksdb
  CONFIGURE_COMMAND ""
  BUILD_COMMAND PORTABLE=1;CFLAGS=-I${CMAKE_CURRENT_SOURCE_DIR}/snappy $(MAKE) -C ${CMAKE_CURRENT_SOURCE_DIR}/rocksdb DISABLE_JEMALLOC=1 static_lib
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(project_rocksdb install_dir)
add_library(rocksdb STATIC IMPORTED)
include_directories(SYSTEM ${install_dir} ${install_dir}/include)
set_property(TARGET rocksdb PROPERTY IMPORTED_LOCATION ${install_dir}/librocksdb.a)
add_dependencies(project_rocksdb project_snappy snappy)
add_dependencies(rocksdb project_rocksdb project_snappy snappy)

# htslib
ExternalProject_Add(
  project_htslib
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/htslib
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/htslib
  CONFIGURE_COMMAND ""
  BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_SOURCE_DIR}/htslib "CFLAGS=-g -Wall -O3" libhts.a
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(project_htslib install_dir)
add_library(htslib STATIC IMPORTED)
include_directories(SYSTEM ${install_dir})
include_directories(SYSTEM ${install_dir}/htslib)

set(htslib_location ${install_dir}/libhts.a)
message(STATUS "htslib target location is ${htslib_location}")
set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${htslib_location})
add_dependencies(htslib project_htslib)


# StatGen
ExternalProject_Add(
  project_statgen
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/statgen
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/statgen
  CONFIGURE_COMMAND ""
  BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_SOURCE_DIR}/statgen USER_WARNINGS=
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(project_statgen install_dir)
add_library(statgen STATIC IMPORTED)
include_directories(SYSTEM ${install_dir})
include_directories(SYSTEM ${install_dir}/include)

set(statgen_location ${install_dir}/libStatGen.a)
message(STATUS "StatGen target location is ${statgen_location}")
set_property(TARGET statgen PROPERTY IMPORTED_LOCATION ${statgen_location})
add_dependencies(statgen project_statgen)

# sparsehash
ExternalProject_Add(
  project_sparsehash
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/sparsehash
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/sparsehash
  CONFIGURE_COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/sparsehash/configure"
  BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_SOURCE_DIR}/sparsehash
  INSTALL_COMMAND ""
)

include_directories(SYSTEM sparsehash/src)

# Find a threads (which will be pthreads on unix)
find_package(Threads)

# ExternalProject_Get_Property(htslib source_dir)
# set(htslib_INCLUDE_DIRS ${source_dir} CACHE INTERNAL "Path to include folder for htslib.")
# set(HTSlib_LIBRARY ${source_dir}/libhts.a CACHE INTERNAL "Path to htslib library files.")

# stations
include_directories(SYSTEM stations/include)

# args
include_directories(SYSTEM args)

# CMake SeqAn
# add_subdirectory(seqan)
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/seqan/include)

# message("SeqAn include path is: ${SEQAN_INCLUDE_PATH}")
# include_directories(SYSTEM ${SEQAN_INCLUDE_DIRS})
# include_directories(SYSTEM ${SEQAN_INCLUDE_PATH})
# set(LIBS ${LIBS} ${SEQAN_LIBRARIES})

# Executables go to a "bin" folder.
set(EXECUTABLE_OUTPUT_PATH ../bin)

# Include public header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# List of all library which require linking
set(graphtyper_all_libraries
  rocksdb
  snappy
  htslib
  statgen
  ${Boost_LIBRARIES}
  pthread
  ${ZLIB_LIBRARIES}
  ${BZIP2_LIBRARIES}
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  message(STATUS "Using GCC")
  list(APPEND graphtyper_all_libraries "rt")
endif ()

# If lz4 is optional for rocksdb
find_package (Lz4)

if (LZ4_FOUND)
  list(APPEND graphtyper_all_libraries ${LZ4_LIBRARY})
endif ()

# If lzma is optional for htslib
find_package (LibLZMA)

if (LIBLZMA_FOUND)
  list(APPEND graphtyper_all_libraries ${LIBLZMA_LIBRARIES})
endif ()

message(STATUS "Libraries: ${graphtyper_all_libraries}")

message(STATUS "Compiling graphtyper's source files")
add_subdirectory(src)

## Testing
enable_testing(true)
add_subdirectory(test)

# add a target to generate API documentation with Doxygen
# find_package(Doxygen)
# if(DOXYGEN_FOUND)
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile @ONLY)
# add_custom_target(docs
#   ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile
#   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
#   COMMENT "Generating API documentation with Doxygen" VERBATIM
# )
# endif(DOXYGEN_FOUND)
