cmake_minimum_required (VERSION 3.1)
project(KMC)

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

find_library(BZIP2 bz2 PATHS ${CMAKE_CURRENT_SOURCE_DIR}/kmc_tools/libs/)
find_library(ZLIB z PATHS ${CMAKE_CURRENT_SOURCE_DIR}/kmc_tools/libs/)

file(GLOB KMER_COUNTER_SOURCES
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/kmer_counter.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/mmer.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/mem_disk_file.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/rev_byte.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/bkb_writer.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/cpu_info.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/bkb_reader.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/fastq_reader.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/timer.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/develop.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/kb_completer.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/kb_storer.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/kmer.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/splitter.cpp
${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/kb_collector.cpp)

add_library(raduls_sse2 STATIC ${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/raduls_sse2.cpp)
target_compile_options(raduls_sse2 PRIVATE "-msse2")
add_library(raduls_sse41 STATIC ${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/raduls_sse41.cpp)
target_compile_options(raduls_sse41 PRIVATE "-msse4.1")
add_library(raduls_avx STATIC ${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/raduls_avx.cpp)
target_compile_options(raduls_avx PRIVATE "-mavx")
add_library(raduls_avx2 STATIC ${CMAKE_CURRENT_SOURCE_DIR}/kmer_counter/raduls_avx2.cpp)
target_compile_options(raduls_avx2 PRIVATE "-mavx2")

file(GLOB KMC_API_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/kmc_api/*.cpp)
add_library(kmcapi STATIC ${KMC_API_SOURCES})

add_library(kmc STATIC ${KMER_COUNTER_SOURCES})

if(NOT APPLE)
  target_link_libraries(kmc kmcapi raduls_sse2 raduls_sse41 raduls_avx raduls_avx2 ${BZIP2} ${ZLIB} -static)
else()
  target_link_libraries(kmc kmcapi raduls_sse2 raduls_sse41 raduls_avx raduls_avx2 ${BZIP2} ${ZLIB})
endif()
