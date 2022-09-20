cmake_minimum_required(VERSION 3.1)
project(KMC VERSION 1.0
  DESCRIPTION "SBWT fork of KMC"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(KMC_MAIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/kmc_core")
set(KMC_CLI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/kmc_CLI")
set(KMC_API_DIR "${CMAKE_CURRENT_SOURCE_DIR}/kmc_api")
set(KMC_TOOLS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/kmc_tools")

find_library(KMC_ZLIB z PATHS ${CMAKE_SOURCE_DIR}/KMC/kmc_tools/libs REQUIRED)
find_library(KMC_BZIP2 bz2 PATHS ${CMAKE_SOURCE_DIR}/KMC/kmc_tools/libs REQUIRED)

if(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -m64 -static-libgcc -static-libstdc++ -pthread -std=c++14")
  set(STATIC_LIBRARY_OPTIONS "-lm -static-libgcc -static-libstdc++ -O3 -pthread -std=c++14")
  add_library(RADULS_SSE2_OBJS OBJECT "${KMC_MAIN_DIR}/raduls_sse2.cpp")
  target_compile_options(RADULS_SSE2_OBJS PRIVATE -msse2)
  add_library(RADULS_SSE41_OBJS OBJECT "${KMC_MAIN_DIR}/raduls_sse41.cpp")
  target_compile_options(RADULS_SSE41_OBJS PRIVATE -msse4.1)
  add_library(RADULS_AVX_OBJS OBJECT "${KMC_MAIN_DIR}/raduls_avx.cpp")
  target_compile_options(RADULS_AVX_OBJS PRIVATE -mavx)
  add_library(RADULS_AVX2_OBJS OBJECT "${KMC_MAIN_DIR}/raduls_avx2.cpp")
  target_compile_options(RADULS_AVX2_OBJS PRIVATE -mavx2)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -m64 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++14")
  set(STATIC_LIBRARY_OPTIONS "-lm -static -O3 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++14")
  add_library(RADULS_SSE2_OBJS OBJECT "${KMC_MAIN_DIR}/raduls_sse2.cpp")
  target_compile_options(RADULS_SSE2_OBJS PRIVATE -msse2)
  add_library(RADULS_SSE41_OBJS OBJECT "${KMC_MAIN_DIR}/raduls_sse41.cpp")
  target_compile_options(RADULS_SSE41_OBJS PRIVATE -msse4.1)
  add_library(RADULS_AVX_OBJS OBJECT "${KMC_MAIN_DIR}/raduls_avx.cpp")
  target_compile_options(RADULS_AVX_OBJS PRIVATE -mavx)
  add_library(RADULS_AVX2_OBJS OBJECT "${KMC_MAIN_DIR}/raduls_avx2.cpp")
  target_compile_options(RADULS_AVX2_OBJS PRIVATE -mavx2)
endif()

add_library(KMC_CLI_OBJS OBJECT "${KMC_CLI_DIR}/kmc.cpp")
add_library(KFF_OBJS OBJECT "${KMC_MAIN_DIR}/kff_writer.cpp")
add_library(KMC_API_OBJS OBJECT "${KMC_API_DIR}/mmer.cpp;${KMC_API_DIR}/kmc_file.cpp;${KMC_API_DIR}/kmer_api.cpp")

add_library(KMC_CORE_OBJS OBJECT "${KMC_MAIN_DIR}/mem_disk_file.cpp;${KMC_MAIN_DIR}/rev_byte.cpp;${KMC_MAIN_DIR}/bkb_writer.cpp;${KMC_MAIN_DIR}/cpu_info.cpp;${KMC_MAIN_DIR}/bkb_reader.cpp;${KMC_MAIN_DIR}/fastq_reader.cpp;${KMC_MAIN_DIR}/timer.cpp;${KMC_MAIN_DIR}/develop.cpp;${KMC_MAIN_DIR}/kb_completer.cpp;${KMC_MAIN_DIR}/kb_storer.cpp;${KMC_MAIN_DIR}/kmer.cpp;${KMC_MAIN_DIR}/splitter.cpp;${KMC_MAIN_DIR}/kb_collector.cpp;${KMC_MAIN_DIR}/kmc_runner.cpp")

add_library(KMC_TOOLS_OBJS OBJECT "${KMC_TOOLS_DIR}/kmer_file_header.cpp;${KMC_TOOLS_DIR}/kmc_tools.cpp;${KMC_TOOLS_DIR}/nc_utils.cpp;${KMC_TOOLS_DIR}/parameters_parser.cpp;${KMC_TOOLS_DIR}/parser.cpp;${KMC_TOOLS_DIR}/tokenizer.cpp;${KMC_TOOLS_DIR}/fastq_filter.cpp;${KMC_TOOLS_DIR}/fastq_reader.cpp;${KMC_TOOLS_DIR}/fastq_writer.cpp;${KMC_TOOLS_DIR}/percent_progress.cpp;${KMC_TOOLS_DIR}/kff_info_reader.cpp;${KMC_API_DIR}/mmer.cpp;${KMC_API_DIR}/kmc_file.cpp;${KMC_API_DIR}/kmer_api.cpp")

add_library(kmc_core STATIC $<TARGET_OBJECTS:KMC_CORE_OBJS> $<TARGET_OBJECTS:KMC_API_OBJS> $<TARGET_OBJECTS:KFF_OBJS> $<TARGET_OBJECTS:RADULS_SSE2_OBJS> $<TARGET_OBJECTS:RADULS_SSE41_OBJS> $<TARGET_OBJECTS:RADULS_AVX_OBJS> $<TARGET_OBJECTS:RADULS_AVX2_OBJS>)

add_library(kmc_tools STATIC $<TARGET_OBJECTS:KMC_TOOLS_OBJS> $<TARGET_OBJECTS:KMC_API_OBJS> $<TARGET_OBJECTS:KFF_OBJS>)
target_link_libraries(kmc_tools PRIVATE ${KMC_ZLIB} ${KMC_BZIP2})

target_include_directories(kmc_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(kmc_tools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
