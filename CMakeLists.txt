cmake_minimum_required (VERSION 3.2)
project (graphtyper)

include(ExternalProject)

# The version number
set (graphtyper_VERSION_MAJOR 2)
set (graphtyper_VERSION_MINOR 7)
set (graphtyper_VERSION_PATCH 2)
set(STATIC_DIR "" CACHE STRING "If set, GraphTyper will be built as a static binary using libraries from the given STATIC_DIR.")

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

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wfatal-errors -pedantic -Wno-variadic-macros -std=c++11 -DSEQAN_HAS_ZLIB=1 -DSEQAN_USE_HTSLIB=1 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0")

if(CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Building in debug mode.")
  set (CMAKE_CXX_FLAGS "-g -O0 -DDEBUG ${CMAKE_CXX_FLAGS}")
elseif(CMAKE_BUILD_TYPE STREQUAL "RELDEB" OR CMAKE_BUILD_TYPE STREQUAL "Reldeb")
  message(STATUS "Building in release (with assertions and debug info) mode.")
  set (CMAKE_CXX_FLAGS "-g -O3 -DDEBUG ${CMAKE_CXX_FLAGS}")
elseif(CMAKE_BUILD_TYPE STREQUAL "PROFILING" OR CMAKE_BUILD_TYPE STREQUAL "Profiling")
  message(STATUS "Building in profiling mode, which is release mode with profiling mode enabled.")
  set (CMAKE_CXX_FLAGS "-pg -O3 -DNDEBUG ${CMAKE_CXX_FLAGS}")
else()
  message(STATUS "Building in release mode.")
  set (CMAKE_CXX_FLAGS "-O3 -DNDEBUG ${CMAKE_CXX_FLAGS}")
endif()

set (CMAKE_CXX_FLAGS_RELEASE "")
set (CMAKE_CXX_FLAGS_DEBUG "")
message(STATUS "CXX flags are: ${CMAKE_CXX_FLAGS}")

# configure a header file to pass some of the CMake settings to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/graphtyper/constants.hpp.in"
  "${PROJECT_BINARY_DIR}/include/graphtyper/constants.hpp"
)

include_directories("${PROJECT_BINARY_DIR}/include")

if (DEFINED CFLAGS)
  set(MYCFLAGS "${CFLAGS}")
  message (STATUS "Custom CFLAGS: ${MYCFLAGS}")
endif()

if (DEFINED LDFLAGS)
  set(MYLDFLAGS "${LDFLAGS}")
  message (STATUS "Custom LDFLAGS: ${MYLDFLAGS}")
endif()

# Find external libraries
## boost
if (NOT DEFINED Boost_USE_STATIC_LIBS)
  set(Boost_USE_STATIC_LIBS ON)
endif()

if (NOT Boost_USE_STATIC_LIBS)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_LOG_DYN_LINK")
endif()

message (STATUS "Checking for Boost: USE_STATIC_LIBS=${Boost_USE_STATIC_LIBS}")
find_package(Boost 1.57.0 REQUIRED COMPONENTS system iostreams log_setup log thread serialization filesystem REQUIRED)
message (STATUS "Boost include directory: ${Boost_INCLUDE_DIR}")
include_directories(SYSTEM ${Boost_INCLUDE_DIR})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ${Boost_INCLUDE_DIR}")

## libdeflate
ExternalProject_Add(
  project_libdeflate
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libdeflate
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/libdeflate
  CONFIGURE_COMMAND ""
  BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_SOURCE_DIR}/libdeflate "CFLAGS=-fPIC -O3 ${CMAKE_C_FLAGS}" libdeflate.a
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(project_libdeflate install_dir)
add_library(libdeflate STATIC IMPORTED)

#include_directories(SYSTEM ${install_dir})

set(libdeflate_location ${install_dir}/libdeflate.a)
message(STATUS "libdeflate target location is ${libdeflate_location}")
set_property(TARGET libdeflate PROPERTY IMPORTED_LOCATION ${libdeflate_location})
add_dependencies(libdeflate project_libdeflate)

## htslib
ExternalProject_Add(
  project_htslib
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/htslib
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/htslib
  CONFIGURE_COMMAND autoheader COMMAND autoconf COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/htslib/configure --disable-libcurl --disable-gcs --with-libdeflate
  "CFLAGS=${MYCFLAGS} -I${CMAKE_CURRENT_SOURCE_DIR}/libdeflate"
  "LDFLAGS=${MYLDFLAGS} -L${CMAKE_CURRENT_SOURCE_DIR}/libdeflate"
  BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_SOURCE_DIR}/htslib "CFLAGS=-g -Wall -O3 ${CMAKE_C_FLAGS} -I${CMAKE_CURRENT_SOURCE_DIR}/libdeflate" libhts.a
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(project_htslib install_dir)
add_library(htslib STATIC IMPORTED)

# both are needed
include_directories(SYSTEM ${install_dir})
include_directories(SYSTEM ${install_dir}/htslib)

set(htslib_location ${install_dir}/libhts.a)
message(STATUS "htslib target location is ${htslib_location}")
set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${htslib_location})
add_dependencies(htslib project_htslib)
add_dependencies(project_htslib libdeflate)

## paw
ExternalProject_Add(
  project_paw
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/paw
  BUILD_IN_SOURCE 1
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/paw
  CONFIGURE_COMMAND ""
  BUILD_COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/paw/build COMMAND ${CMAKE_COMMAND} -H${CMAKE_CURRENT_SOURCE_DIR}/paw -B${CMAKE_CURRENT_SOURCE_DIR}/paw/build -DFORCE_AVX512=1 COMMAND $(MAKE) -C ${CMAKE_CURRENT_SOURCE_DIR}/paw/build static
  INSTALL_COMMAND "")


add_library(paw STATIC IMPORTED)
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/paw/include)
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/paw/build/include)
set(paw_location ${CMAKE_CURRENT_SOURCE_DIR}/paw/build/lib/libpaw.a)
message(STATUS "paw target location is ${paw_location}")
set_property(TARGET paw PROPERTY IMPORTED_LOCATION ${paw_location})
add_dependencies(paw project_paw)

## parallel_hashmap
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/parallel-hashmap)

## Find a multi-thread library (which will likely be pthread on unix)
find_package(Threads)

## SeqAn
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/seqan/include)

## Executables go to a "bin" folder.
set(EXECUTABLE_OUTPUT_PATH ../bin)

## Include public header files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)


if(STATIC_DIR STREQUAL "")
  ### DYNAMIC BUILD
  message(STATUS "Creating dynamic graphtyper binary.")

  ## zlib
  message (STATUS "Checking for zlib")
  find_package(ZLIB REQUIRED)
  include_directories(SYSTEM ${ZLIB_INCLUDE_DIRS})

  ## bzip2
  message (STATUS "Checking for bzip2")
  find_package (BZip2 REQUIRED)
  include_directories(SYSTEM ${BZIP_INCLUDE_DIRS})

  ## List of all library which require linking
  set(graphtyper_all_libraries
    htslib
    paw
    libdeflate
    ${Boost_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${ZLIB_LIBRARIES}
    ${BZIP2_LIBRARIES}
  )

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "Using GCC")
    list(APPEND graphtyper_all_libraries "rt")
  endif ()

  # If lzma is optional for htslib
  find_package (LibLZMA)

  if (LIBLZMA_FOUND)
    list(APPEND graphtyper_all_libraries ${LIBLZMA_LIBRARIES})
  endif ()
else(STATIC_DIR STREQUAL "")
  ### STATIC BUILD
  message(STATUS "Creating static graphtyper binary.")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
  set(graphtyper_all_libraries
    htslib
    paw
    libdeflate
    ${STATIC_DIR}/libboost_system.a
    ${STATIC_DIR}/libboost_iostreams.a
    ${STATIC_DIR}/libboost_log_setup.a
    ${STATIC_DIR}/libboost_log.a
    ${STATIC_DIR}/libboost_thread.a
    ${STATIC_DIR}/libboost_serialization.a
    ${STATIC_DIR}/libboost_filesystem.a
    ${STATIC_DIR}/libboost_regex.a
    ${STATIC_DIR}/libboost_date_time.a
    ${STATIC_DIR}/libboost_chrono.a
    ${STATIC_DIR}/libboost_atomic.a
    ${STATIC_DIR}/libpthread.a
    ${STATIC_DIR}/libz.a
    ${STATIC_DIR}/libbz2.a
    ${STATIC_DIR}/liblzma.a
    ${STATIC_DIR}/librt.a
)
endif(STATIC_DIR STREQUAL "")

message(STATUS "Libraries: ${graphtyper_all_libraries}")
add_subdirectory(src)

## Testing
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/catch/single_include/)

enable_testing(true)
add_subdirectory(test)
