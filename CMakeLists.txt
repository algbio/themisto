cmake_minimum_required (VERSION 3.1)
project(Themisto)
include(ExternalProject)

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

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

## Require c++14
set(CMAKE_CXX_STANDARD 14)
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

## Don't build the parts that we won't need
set_target_properties(divsufsort PROPERTIES EXCLUDE_FROM_ALL 1)

## Set includes
include_directories(${CMAKE_SOURCE_DIR}/zstr ${CMAKE_SOURCE_DIR}/sdsl-lite/include ${CMAKE_SOURCE_DIR}/sdsl-lite/build/external/libdivsufsort/include ${CMAKE_SOURCE_DIR}/KMC/kmer_counter ${CMAKE_SOURCE_DIR}/KMC ${CMAKE_SOURCE_DIR}/KMC/kmc_api ${CMAKE_BINARY_DIR}/include ${CMAKE_SOURCE_DIR}/libwheeler)

## Libraries & executables to build
add_library(kmc_wrapper KMC_wrapper.cpp)
add_executable(pseudoalign pseudoalign.cpp)
add_executable(build_index build_index.cpp)

## Fix compilation order
add_dependencies(kmc_wrapper kmc)
add_dependencies(sdsl divsufsort64)
add_dependencies(pseudoalign kmc_wrapper)
add_dependencies(build_index kmc_wrapper)

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
target_link_libraries(pseudoalign kmc_wrapper sdsl Threads::Threads OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB})
target_link_libraries(build_index kmc_wrapper sdsl Threads::Threads OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB})

## Build tests if doing a Debug build
if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_executable(themisto_tests tests.cpp)
  add_dependencies(themisto_tests kmc_wrapper)
  add_dependencies(sdsl divsufsort64)
  target_link_libraries(themisto_tests kmc_wrapper sdsl Threads::Threads OpenMP::OpenMP_CXX ${BZIP2} ${ZLIB})
else()
  set_target_properties(gtest PROPERTIES EXCLUDE_FROM_ALL 1)
  set_target_properties(gtest_main PROPERTIES EXCLUDE_FROM_ALL 1)
endif()
