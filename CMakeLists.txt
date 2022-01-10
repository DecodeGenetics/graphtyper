cmake_minimum_required (VERSION 3.2)
project (graphtyper)

include(ExternalProject)
include(CheckIncludeFileCXX)

# The version number
set (graphtyper_VERSION_MAJOR 2)
set (graphtyper_VERSION_MINOR 7)
set (graphtyper_VERSION_PATCH 3)

# Graphtyper's headers
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")

## Executables go to a "bin" folder.
set(EXECUTABLE_OUTPUT_PATH ../bin)

#############################################################################
## Build flags
#############################################################################

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wfatal-errors -pedantic -Wno-variadic-macros -std=c++17")

if(CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Building in debug mode.")
  set (CMAKE_CXX_FLAGS "-g -O0 -DDEBUG ${CMAKE_CXX_FLAGS}")
elseif(CMAKE_BUILD_TYPE STREQUAL "RELDEB" OR CMAKE_BUILD_TYPE STREQUAL "Reldeb" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
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


if (DEFINED CFLAGS)
  set(MYCFLAGS "${CFLAGS}")
  message (STATUS "Custom CFLAGS: ${MYCFLAGS}")
endif()

if (DEFINED LDFLAGS)
  set(MYLDFLAGS "${LDFLAGS}")
  message (STATUS "Custom LDFLAGS: ${MYLDFLAGS}")
endif()

#############################################################################
## Passing constants to source
#############################################################################

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

# configure a header file to pass some of the CMake settings to the source code
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/graphtyper/constants.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/graphtyper/constants.hpp"
)

#############################################################################
## Dependency setup
#############################################################################

set (STATIC_DIR "" CACHE STRING "Build in 'static' mode and include libraries in this directory.")

#############################################################################
## htslib
#############################################################################

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/htslib/Makefile)
    message(FATAL_ERROR "htslib not found. Run 'git submodule update --init' to fetch all submodules.")
endif()

ExternalProject_Add(
    project_htslib
    BUILD_IN_SOURCE 1
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/htslib
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/htslib
    CONFIGURE_COMMAND cp -a ${CMAKE_CURRENT_SOURCE_DIR}/htslib ${CMAKE_CURRENT_BINARY_DIR}/ COMMAND autoheader COMMAND autoconf COMMAND ${CMAKE_CURRENT_BINARY_DIR}/htslib/configure --disable-libcurl --disable-gcs --with-libdeflate
        "CFLAGS=${MYCFLAGS} -g -Wall -O3 ${CMAKE_C_FLAGS} -I${CMAKE_CURRENT_BINARY_DIR}/libdeflate"
        "LDFLAGS=${MYLDFLAGS} -L${CMAKE_CURRENT_BINARY_DIR}/libdeflate"
        "CC=${CMAKE_C_COMPILER}"
    BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_BINARY_DIR}/htslib libhts.a
    INSTALL_COMMAND ""
)

add_library(htslib STATIC IMPORTED)

# both are needed
include_directories(SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/htslib)
include_directories(SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/htslib/htslib)

set(htslib_location ${CMAKE_CURRENT_BINARY_DIR}/htslib/libhts.a)
message(STATUS "htslib target location is ${htslib_location}")
set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${htslib_location})
add_dependencies(htslib project_htslib)
add_dependencies(project_htslib libdeflate)
list(APPEND graphtyper_all_libraries "${htslib_location}")

#############################################################################
## paw
#############################################################################

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/paw/CMakeLists.txt)
    message(FATAL_ERROR "paw not found. Run 'git submodule update --init' to fetch all submodules.")
endif()

ExternalProject_Add(
    project_paw
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/paw
    BUILD_IN_SOURCE 0
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/paw
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${CMAKE_COMMAND} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -H${CMAKE_CURRENT_SOURCE_DIR}/paw -B${CMAKE_CURRENT_BINARY_DIR}/paw -DFORCE_AVX512=1
        COMMAND $(MAKE) -C ${CMAKE_CURRENT_BINARY_DIR}/paw static
    INSTALL_COMMAND ""
)

add_library(paw STATIC IMPORTED)
include_directories(SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/paw/include)
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/paw/include)
set(paw_location ${CMAKE_CURRENT_BINARY_DIR}/paw/lib/libpaw.a)
message(STATUS "paw target location is ${paw_location}")
set_property(TARGET paw PROPERTY IMPORTED_LOCATION ${paw_location})
add_dependencies(paw project_paw)
list(APPEND graphtyper_all_libraries "${paw_location}")


#############################################################################
## libdeflate
#############################################################################

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/libdeflate/Makefile)
    message(FATAL_ERROR "libdeflate not found. Run 'git submodule update --init' to fetch all submodules.")
endif ()

ExternalProject_Add(
    project_libdeflate
    BUILD_IN_SOURCE 1
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/libdeflate
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/libdeflate
    CONFIGURE_COMMAND cp -a ${CMAKE_CURRENT_SOURCE_DIR}/libdeflate ${CMAKE_CURRENT_BINARY_DIR}/
    BUILD_COMMAND $(MAKE) -C ${CMAKE_CURRENT_BINARY_DIR}/libdeflate
        PREFIX="${CMAKE_CURRENT_BINARY_DIR}/libdeflate"
        "CC=${CMAKE_C_COMPILER}"
        "CFLAGS=-fPIC -O3 ${CMAKE_C_FLAGS}" libdeflate.a
    INSTALL_COMMAND ""
)

add_library(libdeflate STATIC IMPORTED)
set(libdeflate_location ${CMAKE_CURRENT_BINARY_DIR}/libdeflate/libdeflate.a)
message(STATUS "libdeflate target location is ${libdeflate_location}")

set_property(TARGET libdeflate PROPERTY IMPORTED_LOCATION ${libdeflate_location})
add_dependencies(libdeflate project_libdeflate)
list(APPEND graphtyper_all_libraries "${libdeflate_location}")

#############################################################################
## parallel_hashmap
#############################################################################

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/parallel-hashmap/parallel_hashmap/phmap.h)
    message(FATAL_ERROR "ParallelHashmap not found. Run 'git submodule update --init' to fetch all submodules.")
endif ()

include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/parallel-hashmap/)
message(STATUS "Using ParallelHashMap from ${CMAKE_CURRENT_SOURCE_DIR}/parallel-hashmap.")

#############################################################################
## SeqAn
#############################################################################

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/seqan/include)
    message(FATAL_ERROR "SeqAn not found. Run 'git submodule update --init' to fetch all submodules.")
endif ()

include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/seqan/include)
message(STATUS "Using SeqAn from ${CMAKE_CURRENT_SOURCE_DIR}/seqan/include.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_HAS_ZLIB=1 -DSEQAN_USE_HTSLIB=1 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0")

#############################################################################
## Cereal
#############################################################################

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cereal/include/)
    message(FATAL_ERROR "SeqAn not found. Run 'git submodule update --init' to fetch all submodules.")
endif ()

include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/cereal/include/)
message(STATUS "Using Cereal from ${CMAKE_CURRENT_SOURCE_DIR}/cereal/include.")

#############################################################################
## Compression
#############################################################################

# compression libraries are never "local"
if (STATIC_DIR STREQUAL "")
    message (STATUS "Using system's compression libraries.")

    find_package(ZLIB REQUIRED)
    include_directories(SYSTEM ${ZLIB_INCLUDE_DIRS})
    list(APPEND graphtyper_all_libraries ${ZLIB_LIBRARIES})

    find_package (BZip2 REQUIRED)
    include_directories(SYSTEM ${BZIP_INCLUDE_DIRS})
    list(APPEND graphtyper_all_libraries ${BZIP2_LIBRARIES})

    find_package (LibLZMA)
    if (LIBLZMA_FOUND)
        include_directories(SYSTEM ${LIBLZMA_INCLUDE_DIRS})
        list(APPEND graphtyper_all_libraries ${LIBLZMA_LIBRARIES})
    endif ()
else () # PREBUILT
    message (STATUS "Using prebuilt compression libraries.")
    list(APPEND graphtyper_all_libraries "${STATIC_DIR}/libz.a")
    list(APPEND graphtyper_all_libraries "${STATIC_DIR}/libbz2.a")
    list(APPEND graphtyper_all_libraries "${STATIC_DIR}/liblzma.a")
endif ()

#############################################################################
## System
#############################################################################

if (STATIC_DIR STREQUAL "")
    ## Find a multi-thread library (which will likely be pthread on unix)
    find_package(Threads)
    list(APPEND graphtyper_all_libraries ${CMAKE_THREAD_LIBS_INIT})

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        list(APPEND graphtyper_all_libraries "rt")
        list(APPEND graphtyper_all_libraries "stdc++fs")
    endif ()
else () # PREBUILT
    list(APPEND graphtyper_all_libraries "${STATIC_DIR}/libpthread.a")

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        list(APPEND graphtyper_all_libraries "${STATIC_DIR}/librt.a")
        list(APPEND graphtyper_all_libraries "${STATIC_DIR}/libstdc++fs.a")
    endif ()
endif()

# LOCAL binaries have static GCC, PREBUILT are all static
if (STATIC_DIR STREQUAL "")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
endif ()

#############################################################################
##  DEPENDENCY END
#############################################################################

message(STATUS "Libraries: ${graphtyper_all_libraries}")

# Add graphtyper src
add_subdirectory(src)

#############################################################################
##  TESTING
#############################################################################

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/catch/single_include/)

enable_testing(true)
add_subdirectory(test)

#############################################################################
##  CLANG-FORMAT
#############################################################################

find_program(CLANG_FORMAT "clang-format-12")

if (CLANG_FORMAT STREQUAL "" OR CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
    find_program(CLANG_FORMAT "clang-format")
endif()

if (NOT CLANG_FORMAT STREQUAL "" AND NOT CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
    add_custom_target (check_format "find" "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_SOURCE_DIR}/test" "-name" "\"*.hpp\"" "-o" "-name" "\"*.cpp\"" "-exec" "${CLANG_FORMAT}" "--style=file" "--dry-run" "--Werror" "{}" "+" COMMENT "Checking all files format with clang-format.")
endif()
