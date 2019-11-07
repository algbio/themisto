cmake_minimum_required (VERSION 3.1)
project (Themisto)
include(ExternalProject)

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
endif()
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

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

## Check dependencies
find_package(Threads REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(OpenMP REQUIRED)

## Require c++14
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

## Build KMC during the configure step (it has no CMakeLists.txt)
configure_file(CMakeLists.txt.in KMC/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/KMC )
if(result)
  message(FATAL_ERROR "CMake step for KMC failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/KMC )
if(result)
  message(FATAL_ERROR "Build step for KMC failed: ${result}")
endif()

## Add local dependencies as targets
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite
                 ${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite/build)
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/BD_BWT_index
                 ${CMAKE_CURRENT_BINARY_DIR}/BD_BWT_index/build)

## Set includes
include_directories(BD_BWT_index/include sdsl-lite/include KMC KMC/kmer_counter)

## Libraries & executables to build
add_library(kmc_wrapper STATIC KMC_wrapper.cpp)
add_executable(pseudoalign pseudoalign.cpp)
add_executable(build_index build_index.cpp)
add_executable(themisto_tests tests.cpp)

## Fix compilation order
add_dependencies(pseudoalign sdsl bdbwt kmc_wrapper divsufsort64)
add_dependencies(build_index sdsl bdbwt kmc_wrapper divsufsort64)
add_dependencies(themisto_tests sdsl bdbwt kmc_wrapper divsufsort64)
add_dependencies(bdbwt sdsl divsufsort64 gtest gtest_main)
add_dependencies(tests bdbwt sdsl divsufsort64 gtest gtest_main)

## Glob the KMC libraries that we need for linking
file(GLOB KMC_OBJECTS KMC/kmer_counter/mmer.o KMC/kmer_counter/mem_disk_file.o KMC/kmer_counter/rev_byte.o KMC/kmer_counter/bkb_writer.o KMC/kmer_counter/cpu_info.o KMC/kmer_counter/bkb_reader.o KMC/kmer_counter/fastq_reader.o KMC/kmer_counter/timer.o KMC/kmer_counter/develop.o KMC/kmer_counter/kb_completer.o KMC/kmer_counter/kb_storer.o KMC/kmer_counter/kmer.o KMC/kmer_counter/splitter.o KMC/kmer_counter/kb_collector.o KMC/kmer_counter/raduls_sse2.o KMC/kmer_counter/raduls_sse41.o KMC/kmer_counter/raduls_avx2.o KMC/kmer_counter/raduls_avx.o KMC/kmc_api/kmer_api.o KMC/kmc_api/kmc_file.o)

## Tell CMake KMC libs are already built
set(OBJS ${KMC_OBJECTS})
set_source_files_properties(${OBJS} PROPERTIES EXTERNAL_OBJECT true GENERATED true)

find_library(BZIP2 bz2 PATHS ${CMAKE_SOURCE_DIR}/KMC/kmc_tools/libs/)
find_library(ZLIB z PATHS ${CMAKE_SOURCE_DIR}/KMC/kmc_tools/libs/)

## Link
target_link_libraries(kmc_wrapper Threads::Threads OpenMP::OpenMP_CXX ${OBJS} ${BZIP2} ${ZLIB})
target_link_libraries(pseudoalign kmc_wrapper sdsl bdbwt divsufsort64 Threads::Threads OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB})
target_link_libraries(build_index kmc_wrapper sdsl bdbwt divsufsort64 Threads::Threads OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB})
target_link_libraries(themisto_tests kmc_wrapper sdsl bdbwt divsufsort64 OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB})
