cmake_minimum_required (VERSION 3.13)
project(Themisto)
include(ExternalProject)

## Set maximum k-mer length
if(NOT MAX_KMER_LENGTH)
  set(MAX_KMER_LENGTH 32) ## Default max k-mer length is 32
endif()
if(BUILD_THEMISTO_TESTS)
  set(MAX_KMER_LENGTH 255) # 255 required in tests
endif()

if(MAX_KMER_LENGTH GREATER 255)
  message(${MAX_KMER_LENGTH})
  message(FATAL_ERROR "Error: MAX_KMER_LENGTH can be at most 255")
endif()

message(STATUS "Maximum k-mer length is: " ${MAX_KMER_LENGTH})
message(STATUS "To change the maximum k-mer length, use the option -DMAX_KMER_LENGTH=n")
message(STATUS "(n can be up to 255)")

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
endif()

# https://stackoverflow.com/questions/48825416/missing-openmp-c-flags-openmp-c-lib-names
if(APPLE)
if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    set(OpenMP_C "${CMAKE_C_COMPILER}" CACHE STRING "" FORCE)
    set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument" CACHE STRING "" FORCE)
    set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5" CACHE STRING "" FORCE)
    set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES} CACHE STRING "" FORCE)
    set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES} CACHE STRING "" FORCE)
    set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES} CACHE STRING "" FORCE)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(OpenMP_CXX "${CMAKE_CXX_COMPILER}" CACHE STRING "" FORCE)
  set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument" CACHE STRING "" FORCE)
  set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5" CACHE STRING "" FORCE)
  set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES} CACHE STRING "" FORCE)
  set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES} CACHE STRING "" FORCE)
  set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES} CACHE STRING "" FORCE)
endif()
endif()

# -Wno-deprecated-copy is for warnings in Roaring.
# -Wno-deprecated-declarations is for warnings in SDSL.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-sign-compare -Wno-deprecated-copy -Wno-deprecated-declarations")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Our sources, excluding the file with the main function
set(THEMISTO_SOURCES
  src/coloring/Roaring_Color_Set.cpp
  src/coloring/coloring.cpp
  src/coloring/color_set.cpp
  src/coloring/color_set_diagnostics.cpp
  src/dump_color_matrix_main.cpp
  src/pseudoalign_main.cpp
  src/pseudoalign.cpp
  src/globals.cpp
  src/run_tracker.cpp
  src/test_tools.cpp
  src/WorkDispatcher.cpp
  src/zpipe.cpp
  src/build_index_main.cpp
  src/extract_unitigs_main.cpp
  src/new_extract_unitigs.cpp
  src/DBG.cpp
  src/stats_main.cpp
  src/make_d_equal_1.cpp
  src/dump_distinct_color_sets_to_binary.cpp
  )

  ## Require zlib
if(NOT ZLIB)
  if(NOT CMAKE_BUILD_ZLIB)
    find_library(ZLIB NAMES z REQUIRED)
  endif()

  if(CMAKE_BUILD_ZLIB OR ZLIB MATCHES NOTFOUND)
    configure_file(CMakeLists_zlib.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/zlib-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/zlib-download )
    if(result)
      message(FATAL_ERROR "CMake step for zlib failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/zlib-download )
    if(result)
      message(FATAL_ERROR "Build step for zlib failed: ${result}")
    endif()
    set(ZLIB ${CMAKE_CURRENT_BINARY_DIR}/external/zlib/libz.a)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/external/zlib)
    ExternalProject_Add(zlib
      SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/external/zlib"
      BUILD_IN_SOURCE   1
      BUILD_COMMAND     $(MAKE) -f ${CMAKE_CURRENT_BINARY_DIR}/external/zlib/Makefile CFLAGS=${CMAKE_C_FLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}
      CONFIGURE_COMMAND ""
      INSTALL_COMMAND   ""
      TEST_COMMAND      ""
      UPDATE_COMMAND    ""
      )
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/external/zlib)
  endif()
endif()

## Require c++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fix for GCC (see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90050 )
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CXX_FILESYSTEM_LIBRARIES "stdc++fs")
endif()

## Check dependencies
find_package(Threads REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(OpenMP REQUIRED)

## Get version number from git tags
find_package(Git)
if(GIT_FOUND)
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
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
add_subdirectory(${CMAKE_SOURCE_DIR}/SBWT
		 ${CMAKE_BINARY_DIR}/external/SBWT/build)

add_custom_target(
   ggcat_cpp_api
   COMMAND make
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ggcat/crates/capi/ggcat-cpp-api
)

set(LIBSBWT ${CMAKE_BINARY_DIR}/external/SBWT/build/lib/libsbwt_static.a)

include(FetchContent)

FetchContent_Declare(
  roaring
  GIT_REPOSITORY https://github.com/RoaringBitmap/CRoaring.git
  GIT_TAG v0.2.66
  GIT_SHALLOW TRUE)

set(ENABLE_ROARING_TESTS OFF CACHE INTERNAL "")

set(ROARING_BUILD_STATIC ON CACHE INTERNAL "")
FetchContent_MakeAvailable(roaring)

FetchContent_GetProperties(roaring)
SET(CPP_ROARING_HEADERS ${roaring_SOURCE_DIR}/cpp/roaring64map.hh  ${roaring_SOURCE_DIR}/cpp/roaring.hh)
file(COPY  ${CPP_ROARING_HEADERS} DESTINATION ${roaring_SOURCE_DIR}/include/roaring)

## Set includes
include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/include/libwheeler
  ${CMAKE_SOURCE_DIR}/zstr
  ${CMAKE_SOURCE_DIR}/SBWT/sdsl-lite/include
  ${CMAKE_SOURCE_DIR}/build/external/SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include
  ${CMAKE_SOURCE_DIR}/SBWT/KMC/kmer_counter
  ${CMAKE_SOURCE_DIR}/SBWT
  ${CMAKE_SOURCE_DIR}/SBWT/KMC
  ${CMAKE_SOURCE_DIR}/SBWT/include
  ${CMAKE_SOURCE_DIR}/SBWT/include/sbwt
  ${CMAKE_SOURCE_DIR}/SBWT/SeqIO/include
  ${CMAKE_SOURCE_DIR}/SBWT/KMC/kmc_api
  ${CMAKE_BINARY_DIR}/include
  ${CMAKE_SOURCE_DIR}/SBWT/KMC/kmc_core
  ${CMAKE_BINARY_DIR}/_deps/roaring-src/include/
  ${CMAKE_SOURCE_DIR}/ggcat/crates/capi/ggcat-cpp-api/include)

set(GGCAT_CPP_BINDINGS ${CMAKE_SOURCE_DIR}/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_cpp_bindings.a)
set(GGCAT_CXX_INTEROP ${CMAKE_SOURCE_DIR}/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_cxx_interop.a)
set(GGCAT_API ${CMAKE_SOURCE_DIR}/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_api.a)

# Build main executable
add_executable(themisto src/themisto_main.cpp ${THEMISTO_SOURCES})
target_compile_definitions(themisto PUBLIC MAX_KMER_LENGTH=${MAX_KMER_LENGTH}) # Define for compiler.
# Fix compilation order
add_dependencies(sdsl divsufsort64)
add_dependencies(themisto ggcat_cpp_api sbwt_static)
## Link
target_link_libraries(themisto PRIVATE
  sdsl
  Threads::Threads
  OpenMP::OpenMP_CXX
  sbwt_static
  ${GGCAT}
  ${ZLIB}
  ${CXX_FILESYSTEM_LIBRARIES}
  kmc_tools
  kmc_core
  roaring
  ${GGCAT_API}
  ${GGCAT_CPP_BINDINGS}
  ${GGCAT_CXX_INTEROP}
  ${CMAKE_DL_LIBS})

# Build tests if requested
if (BUILD_THEMISTO_TESTS)
  add_subdirectory(${CMAKE_SOURCE_DIR}/googletest
    ${CMAKE_BINARY_DIR}/external/googletest/build)
  add_executable(themisto_tests tests/test_main.cpp ${THEMISTO_SOURCES})
  target_compile_definitions(themisto_tests PUBLIC MAX_KMER_LENGTH=${MAX_KMER_LENGTH}) # Define for compiler.
  add_dependencies(sdsl divsufsort64)
  add_dependencies(themisto_tests ggcat_cpp_api sbwt_static)
  target_include_directories(themisto_tests PRIVATE ${CMAKE_SOURCE_DIR}/googletest/googletest/include)
  target_link_libraries(themisto_tests PRIVATE
    gtest sdsl
    Threads::Threads
    OpenMP::OpenMP_CXX
    sbwt_static
    ${GGCAT}
    ${ZLIB}
    ${CXX_FILESYSTEM_LIBRARIES}
    kmc_tools
    kmc_core
    roaring
    ${GGCAT_API}
    ${GGCAT_CPP_BINDINGS}
    ${GGCAT_CXX_INTEROP}
    ${CMAKE_DL_LIBS})
endif()

if(BUILD_IO_BENCHMARK)
  message("Setting up IO benchmark.")
  add_executable(benchmark_io tests/benchmark_io.cpp ${THEMISTO_SOURCES})
  target_link_libraries(benchmark_io PRIVATE sdsl
    Threads::Threads
    OpenMP::OpenMP_CXX
    ${ZLIB}
    ${CXX_FILESYSTEM_LIBRARIES}
    kmc_tools
    kmc_core
    roaring)
endif()
