# Installation
## Requirements
External dependencies: libz and libbz2.

Other requirements: C++14 compliant compiler, CMake v3.10 or newer.

## Compiling
First, install libz and libbz2 if your system does not already have
those. Then run:

	cmake .
    make

This will produce the build\_index, pseudoalignment, and
themisto\_tests executables.

(Optional) Check if everything installed correctly by running tests:

    ./themisto_tests

If there is a problem, the tests will terminate and report an error message.

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
