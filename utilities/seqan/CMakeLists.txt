cmake_minimum_required(VERSION 2.8.8)
message("Fetching SeqAn from Github")
project(seqan_builder CXX)
include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
    seqan
    PREFIX ${CMAKE_BINARY_DIR}/seqan
    GIT_REPOSITORY https://github.com/seqan/seqan.git
    GIT_TAG develop
    TIMEOUT 30
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull origin develop
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
   )

# Expose required variable (SEQAN_INCLUDE_PATH) to parent scope
ExternalProject_Get_Property(seqan source_dir)
set(SEQAN_INCLUDE_PATH ${source_dir}/include CACHE INTERNAL "Path to include folder for Catch")
