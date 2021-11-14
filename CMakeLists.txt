cmake_minimum_required (VERSION 3.1)
project(Themisto)
include(ExternalProject)

## Set maximum k-mer length
if(NOT MAX_KMER_LENGTH)
  set(MAX_KMER_LENGTH 63) ## Default max k-mer length is 63
endif()
if(MAX_KMER_LENGTH GREATER 255)
  message(${MAX_KMER_LENGTH})
  message(FATAL_ERROR "Error: MAX_KMER_LENGTH can be at most 255")
endif()
MATH(EXPR K_PLUS_1_MER_LENGTH "${MAX_KMER_LENGTH}+1") # +1 because the construction requires (k+1)-mers
add_definitions(-DKMER_MAX_LENGTH=${K_PLUS_1_MER_LENGTH}) # Define for compiler. Note that the variable name here is KMER_MAX_LENGTH, not MAX_KMER_LENGTH, but that is correct.
message(STATUS "Maximum k-mer length is: " ${MAX_KMER_LENGTH})
message(STATUS "To change the maximum k-mer length, use the option -DMAX_KMER_LENGTH=n")
message(STATUS "(n can be up to 255)")

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Our sources, excluding the file with the main function
set(THEMISTO_SOURCES src/pseudoalign_main.cpp src/globals.cpp src/EM_sort.cpp src/Block.cpp src/EM_algorithms.cpp src/test_tools.cpp src/Themisto.cpp src/WorkDispatcher.cpp src/zpipe.cpp src/build_index_main.cpp src/pseudoalign_main.cpp)

## Check if we should compile zlib or bzip2
if(CMAKE_BUILD_ZLIB)
  configure_file(CMakeLists_zlib.txt.in ${CMAKE_BINARY_DIR}/external/zlib-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/zlib-download )
  if(result)
    message(FATAL_ERROR "CMake step for zlib failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/zlib-download )
  if(result)
    message(FATAL_ERROR "Build step for zlib failed: ${result}")
  endif()
  set(ZLIB ${CMAKE_BINARY_DIR}/lib/libz.a)
  include_directories(${CMAKE_BINARY_DIR}/external/zlib)
else()
  find_library(ZLIB z PATHS ${CMAKE_SOURCE_DIR}/KMC/kmc_tools/libs)
endif()
if(CMAKE_BUILD_BZIP2)
  configure_file(CMakeLists_bzip2.txt.in ${CMAKE_BINARY_DIR}/external/bzip2-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/bzip2-download )
  if(result)
    message(FATAL_ERROR "CMake step for bzip2 failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/bzip2-download )
  if(result)
    message(FATAL_ERROR "Build step for bzip2 failed: ${result}")
  endif()
  set(BZIP2 ${CMAKE_BINARY_DIR}/external/bzip2/libbz2.a)
  ExternalProject_Add(bzip2
    SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/external/bzip2"
    BUILD_IN_SOURCE   1
    BUILD_COMMAND     $(MAKE) -f ${CMAKE_CURRENT_BINARY_DIR}/external/bzip2/Makefile CFLAGS=${CMAKE_C_FLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    TEST_COMMAND      ""
    UPDATE_COMMAND    ""
  )
  include_directories(${CMAKE_BINARY_DIR}/external/bzip2)
else()
  find_library(BZIP2 bz2 PATHS ${CMAKE_BINARY_DIR}/external/bzip2)
endif()

## Require c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## Check dependencies
find_package(Threads REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(OpenMP REQUIRED)

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
configure_file(${CMAKE_SOURCE_DIR}/version.h.in ${CMAKE_BINARY_DIR}/include/version.h @ONLY)

## Add local dependencies as targets
add_subdirectory(${CMAKE_SOURCE_DIR}/sdsl-lite
                 ${CMAKE_BINARY_DIR}/external/sdsl-lite/build)
add_subdirectory(${CMAKE_SOURCE_DIR}/KMC
                 ${CMAKE_BINARY_DIR}/external/KMC/build)

## Set includes
include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/libwheeler ${CMAKE_SOURCE_DIR}/zstr ${CMAKE_SOURCE_DIR}/sdsl-lite/include ${CMAKE_SOURCE_DIR}/build/external/sdsl-lite/build/external/libdivsufsort/include ${CMAKE_SOURCE_DIR}/KMC/kmer_counter ${CMAKE_SOURCE_DIR}/KMC ${CMAKE_SOURCE_DIR}/KMC/kmc_api ${CMAKE_BINARY_DIR}/include)

## Libraries & executables to build
add_library(kmc_wrapper src/KMC_wrapper.cpp)
add_executable(themisto src/themisto_main.cpp ${THEMISTO_SOURCES})

## Fix compilation order
add_dependencies(kmc_wrapper kmc)
add_dependencies(sdsl divsufsort64)
add_dependencies(themisto kmc_wrapper)

## Handle compiling zlib or bzip2 from source
if(CMAKE_BUILD_ZLIB)
  add_subdirectory(${CMAKE_BINARY_DIR}/external/zlib
                   ${CMAKE_BINARY_DIR}/external/zlib)
  set_target_properties(example PROPERTIES EXCLUDE_FROM_ALL 1)
  set_target_properties(minigzip PROPERTIES EXCLUDE_FROM_ALL 1)
  add_dependencies(kmc_wrapper zlibstatic)
endif()
if(CMAKE_BUILD_BZIP2)
  add_dependencies(kmc_wrapper bzip2)
endif()

## Link
target_link_libraries(kmc_wrapper kmc Threads::Threads OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB})
target_link_libraries(themisto kmc_wrapper sdsl Threads::Threads OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB})

## Build tests if doing a Debug build
if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_executable(themisto_tests tests/test_main.cpp ${THEMISTO_SOURCES})

  add_dependencies(themisto_tests kmc_wrapper)

  add_dependencies(sdsl divsufsort64)
  find_library(GTEST NAMES gtest PATHS googletest/build/lib/)
  target_link_libraries(themisto_tests kmc_wrapper ${GTEST} sdsl Threads::Threads OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB})
else()
#  set_target_properties(gtest PROPERTIES EXCLUDE_FROM_ALL 1)
#  set_target_properties(gtest_main PROPERTIES EXCLUDE_FROM_ALL 1)
endif()
