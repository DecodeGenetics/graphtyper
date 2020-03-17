cmake_minimum_required (VERSION 2.8.11)
project (graphtyper)

include(ExternalProject)

# The version number.
set (graphtyper_VERSION_MAJOR 2)
set (graphtyper_VERSION_MINOR 2)
set (graphtyper_VERSION_PATCH 1)

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated SHA commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_SHORT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest SHA commit hash of the working branch
execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_LONG_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest SHA commit hash of the working branch
execute_process(
  COMMAND git diff COMMAND wc -l
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_NUM_DIRTY_LINES
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if ($ENV{GT_BUILD_TYPE})
  set(GT_BUILD_TYPE $ENV{GT_BUILD_TYPE})
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wfatal-errors -pedantic -Wno-variadic-macros -std=c++11 -DSEQAN_HAS_ZLIB=1 -DSEQAN_USE_HTSLIB=1 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0")

if(CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Building in debug mode.")
  set (CMAKE_CXX_FLAGS "-g -O0 -DDEBUG ${CMAKE_CXX_FLAGS}")
elseif(GT_BUILD_TYPE STREQUAL "RELDEB" OR GT_BUILD_TYPE STREQUAL "Reldeb")
  message(STATUS "Building in release (with assertions and debug info) mode.")
  set (CMAKE_CXX_FLAGS "-g -O3 -DDEBUG ${CMAKE_CXX_FLAGS}")
elseif(GT_BUILD_TYPE STREQUAL "PROFILING" OR GT_BUILD_TYPE STREQUAL "Profiling")
  message(STATUS "Building in profiling mode, which is release mode with profiling mode enabled.")
  set (CMAKE_CXX_FLAGS "-pg -O3 -DNDEBUG ${CMAKE_CXX_FLAGS}")
else()
  message(STATUS "Building in release mode.")
  set (CMAKE_CXX_FLAGS "-O3 -DNDEBUG ${CMAKE_CXX_FLAGS}")
endif()

message(STATUS "CXX flags are: ${CMAKE_CXX_FLAGS}")

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/utilities/cmake_modules/")

# configure a header file to pass some of the CMake settings to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/graphtyper/constants.hpp.in"
  "${PROJECT_BINARY_DIR}/include/graphtyper/constants.hpp"
)

include_directories("${PROJECT_BINARY_DIR}/include")

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
find_package(Boost 1.57.0 REQUIRED COMPONENTS system iostreams log_setup log thread serialization filesystem REQUIRED)
message (STATUS "Boost include directory: ${Boost_INCLUDE_DIR}")
include_directories(SYSTEM ${Boost_INCLUDE_DIR})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${Boost_INCLUDE_DIR}")

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
  CONFIGURE_COMMAND autoheader COMMAND autoconf COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/htslib/configure --disable-libcurl
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

# Find a multi-thread library (which will likely be pthread on unix)
find_package(Threads)

# ExternalProject_Get_Property(htslib source_dir)
# set(htslib_INCLUDE_DIRS ${source_dir} CACHE INTERNAL "Path to include folder for htslib.")
# set(HTSlib_LIBRARY ${source_dir}/libhts.a CACHE INTERNAL "Path to htslib library files.")

# paw
ExternalProject_Add(
  project_paw
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/paw
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/paw
  CONFIGURE_COMMAND ""
  BUILD_COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/paw/build COMMAND cmake -H${CMAKE_CURRENT_SOURCE_DIR}/paw -B${CMAKE_CURRENT_SOURCE_DIR}/paw/build COMMAND $(MAKE) -C ${CMAKE_CURRENT_SOURCE_DIR}/paw/build static
  INSTALL_COMMAND ""
)

add_library(paw STATIC IMPORTED)
include_directories(SYSTEM paw/include)
set(paw_location ${CMAKE_CURRENT_SOURCE_DIR}/paw/build/lib/libpaw.a)
message(STATUS "paw target location is ${paw_location}")
set_property(TARGET paw PROPERTY IMPORTED_LOCATION ${paw_location})
add_dependencies(paw project_paw)

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
  paw
  ${Boost_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${ZLIB_LIBRARIES}
  ${BZIP2_LIBRARIES}
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  message(STATUS "Using GCC")
  list(APPEND graphtyper_all_libraries "rt")
endif ()

# lz4 is optional for rocksdb
find_package (Lz4)

if (LZ4_FOUND)
  list(APPEND graphtyper_all_libraries ${LZ4_LIBRARY})
endif ()

# If lzma is optional for htslib
find_package (LibLZMA)

if (LIBLZMA_FOUND)
  list(APPEND graphtyper_all_libraries ${LIBLZMA_LIBRARIES})
endif ()

# zstd is optional for rocksdb
find_package (ZStd)

if (ZSTD_FOUND)
  list(APPEND graphtyper_all_libraries ${ZSTD_LIBRARIES})
endif()


message(STATUS "Libraries: ${graphtyper_all_libraries}")
add_subdirectory(src)

## Testing
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/catch/single_include/)

enable_testing(true)
add_subdirectory(test)
