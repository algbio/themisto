External dependencies: libz and libbz2

Compiling (todo: better instructions and/or build process):

First, install libz and libbz2 if your system does not already have those. Then:

    # Compile sdsl-lite
    cd sdsl-lite
    sh install.sh
    
    # Compile BD_BWT_index
    cd ../BD_BWT_index
    cmake .
    make

    # Compile KMC
    cd ../KMC
    make
    cd ..

    # Compile the main executables
    make tests build_index pseudoalign

To run tests, run:

    ./tests

If there is a problem, it will terminate and report an error message.

To build the index, run:

    ./build_index

The program will then print instructions on how to use it.

To pseudoalign against the index, run:

    ./pseudoalign

The program will then print instructions on how to use it.
