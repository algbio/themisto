cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
project(BD_BWT_index)

file(GLOB LIBRARY_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/dbwt/*.c src/io_tools.cpp)
file(GLOB TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

include_directories(include dbwt ${CMAKE_SOURCE_DIR}/sdsl-lite/include ${CMAKE_SOURCE_DIR}/sdsl-lite/build/external/libdivsufsort/include)

add_library(bdbwt ${LIBRARY_SOURCES} ${COMMON_SOURCES})
target_link_libraries(bdbwt sdsl divsufsort64)

add_executable(bdbwt_tests ${TEST_SOURCES} ${COMMON_SOURCES})
target_link_libraries(bdbwt_tests sdsl divsufsort64)
