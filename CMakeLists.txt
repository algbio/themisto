cmake_minimum_required (VERSION 3.1)
project(Themisto)

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

## Require c++14
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

## Check dependencies
find_package(Threads REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(OpenMP REQUIRED)

find_library(BZIP2 bz2 HINTS ${CMAKE_SOURCE_DIR}/KMC/kmc_tools/libs)
find_library(ZLIB z HINTS ${CMAKE_SOURCE_DIR}/KMC/kmc_tools/libs)

## Get version number from git tags
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags
    WORKING_DIRECTORY "${local_dir}"
    OUTPUT_VARIABLE THEMISTO_BUILD_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message( STATUS "Themisto version: ${THEMISTO_BUILD_VERSION}")
else()
  message(STATUS "GIT not found, version numbering disabled")
  set(THEMISTO_BUILD_VERSION "unknown")
endif()
string(TIMESTAMP THEMISTO_BUILD_TIMESTAMP)

## Generate a version.h file containing build version and timestamp
configure_file(${CMAKE_SOURCE_DIR}/version.h.in ${CMAKE_SOURCE_DIR}/version.h @ONLY)

## Add local dependencies as targets
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite
                 ${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite/build)
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/BD_BWT_index
                 ${CMAKE_CURRENT_BINARY_DIR}/BD_BWT_index/build)
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/KMC
                 ${CMAKE_CURRENT_BINARY_DIR}/KMC/build)

## Set includes
include_directories(BD_BWT_index/include sdsl-lite/include KMC/kmer_counter KMC)

## Libraries & executables to build
add_library(kmc_wrapper STATIC KMC_wrapper.cpp)
add_executable(pseudoalign pseudoalign.cpp)
add_executable(build_index build_index.cpp)
add_executable(themisto_tests tests.cpp)

## Fix compilation order
add_dependencies(kmc_wrapper kmc)
add_dependencies(bdbwt sdsl divsufsort64)
add_dependencies(bdbwt_tests sdsl divsufsort64)
add_dependencies(pseudoalign kmc_wrapper bdbwt)
add_dependencies(build_index kmc_wrapper bdbwt)
add_dependencies(themisto_tests kmc_wrapper bdbwt)

## Link
target_link_libraries(kmc_wrapper kmc Threads::Threads OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB} -static)
target_link_libraries(pseudoalign kmc_wrapper bdbwt Threads::Threads OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB})
target_link_libraries(build_index kmc_wrapper bdbwt Threads::Threads OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB})
target_link_libraries(themisto_tests kmc_wrapper bdbwt Threads::Threads OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB})
