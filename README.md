# Installation
## Requirements
Compilation: C++14 compliant compiler with OpenMP support, and CMake v3.1 or newer.

## Compiling
Enter the Themisto directory and run

	cd build
	cmake ..
    make

This will produce the build\_index, pseudoalignment, and
themisto\_tests executables in the build/bin/ directory.

If you run into problems involving zlib or bzip2, you can instruct the
build process to download & compile them from source with

	cmake -DCMAKE_BUILD_ZLIB=1 -DCMAKE_BUILD_BZIP2=1 ..
	make

## Compiling on macOS
Compiling Themisto on macOS requires users to first install gcc-6.1 or
newer from homebrew with

	brew install gcc@6

Afterwards, Themisto can be compiled by entering the directory and running

	cd build
	cmake -DCMAKE_C_COMPILER=$(which gcc-6) -DCMAKE_CXX_COMPILER=$(which g++-6) ..
	make

Note that macOS has a very small limit for the number of concurrently
opened files. Themisto can use temporary files to conserve RAM, and
may run into this limit. To increase the limit, run the command

	ulimit -n 2048

# Usage
## Indexing
To build the index, run:

    build_index

The program will then print instructions on how to use it.

## Pseudoalignment
To pseudoalign against the index, run:

    pseudoalign

The program will then print instructions on how to use it.

# License

This software is licensed under GPLv2. See LICENSE.txt.
