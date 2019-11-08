cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

project(BD_BWT_index)

set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

file(GLOB LIBRARY_SOURCES dbwt/*.c src/io_tools.cpp)
file(GLOB TEST_SOURCES src/main.cpp)

include_directories(include dbwt ${CMAKE_SOURCE_DIR}/sdsl-lite/include ${CMAKE_SOURCE_DIR}/sdsl-lite/build/external/libdivsufsort/include)

add_library(bdbwt ${LIBRARY_SOURCES} ${COMMON_SOURCES})
target_link_libraries(bdbwt sdsl divsufsort64)

add_executable(tests ${TEST_SOURCES} ${COMMON_SOURCES})
target_link_libraries(tests sdsl divsufsort64)
