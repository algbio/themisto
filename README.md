# Installation
## Requirements
Compilation: C++14 compliant compiler with OpenMPI support, and CMake v3.1 or newer.

## Compiling
Enter the Themisto directory and run

	cmake .
    make

This will produce the build\_index, pseudoalignment, and
themisto\_tests executables.

(Optional) Check if everything installed correctly by running tests:

    ./themisto_tests

If there is a problem, the tests will terminate and report an error message.

## Compiling on macOS
Compiling Themisto on macOS requires users to first install gcc-6.1 or
newer from homebrew with

	brew install gcc@6

Afterwards, Themisto can be compiled by entering the directory and running

	cmake -DCMAKE_C_COMPILER=$(which gcc-6) -DCMAKE_CXX_COMPILER=$(which g++-6) .
	make

# Usage
## Indexing
To build the index, run:

    ./build_index

The program will then print instructions on how to use it.

## Pseudoalignment
To pseudoalign against the index, run:

    ./pseudoalign

The program will then print instructions on how to use it.

# License

This software is licensed under GPLv2. See LICENSE.txt.
