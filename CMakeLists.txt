cmake_minimum_required (VERSION 3.1)
project (Themisto)
include(ExternalProject)

## Check dependencies
find_package(Threads REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

## Add local dependencies as targets
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite
                 ${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite/build)
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/BD_BWT_index
                 ${CMAKE_CURRENT_BINARY_DIR}/BD_BWT_index/build)
ExternalProject_Add(KMC
  SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/KMC"
  BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/KMC"
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND   ""
  UPDATE_COMMAND    ""
  TEST_COMMAND      ""
)

## Set includes
include_directories(BD_BWT_index/include sdsl-lite/include KMC KMC/kmer_counter/)

## Glob all libraries that we need for linking
file(GLOB OWN_LIBRARIES BD_BWT_index/lib/*.a sdsl-lite/build/lib/libsdsl.a sdsl-lite/build/external/libdivsufsort/lib/libdivsufsort64.a)
file(GLOB KMC_OBJECTS KMC/kmer_counter/mmer.o KMC/kmer_counter/mem_disk_file.o KMC/kmer_counter/rev_byte.o KMC/kmer_counter/bkb_writer.o KMC/kmer_counter/cpu_info.o KMC/kmer_counter/bkb_reader.o KMC/kmer_counter/fastq_reader.o KMC/kmer_counter/timer.o KMC/kmer_counter/develop.o KMC/kmer_counter/kb_completer.o KMC/kmer_counter/kb_storer.o KMC/kmer_counter/kmer.o KMC/kmer_counter/splitter.o KMC/kmer_counter/kb_collector.o KMC/kmer_counter/raduls_sse2.o KMC/kmer_counter/raduls_sse41.o KMC/kmer_counter/raduls_avx2.o KMC/kmer_counter/raduls_avx.o KMC/kmc_api/kmer_api.o KMC/kmc_api/kmc_file.o)

add_library(kmc_wrapper STATIC KMC_wrapper.cpp)
target_link_libraries(kmc_wrapper ZLIB::ZLIB bz2 Threads::Threads OpenMP::OpenMP_CXX)

add_executable(pseudoalign pseudoalign.cpp)
target_link_libraries(pseudoalign kmc_wrapper ZLIB::ZLIB bz2 Threads::Threads OpenMP::OpenMP_CXX ${OWN_LIBRARIES} ${KMC_OBJECTS})

add_executable(build_index build_index.cpp)
target_link_libraries(build_index kmc_wrapper ZLIB::ZLIB bz2 Threads::Threads OpenMP::OpenMP_CXX ${OWN_LIBRARIES} ${KMC_OBJECTS})

add_executable(themisto_tests tests.cpp)
target_link_libraries(themisto_tests kmc_wrapper ZLIB::ZLIB bz2 OpenMP::OpenMP_CXX ${OWN_LIBRARIES} ${KMC_OBJECTS})
