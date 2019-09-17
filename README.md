KallistoLite is working title.

Compiling:

    cd sdsl-lite
    sh install.sh
    cd ..
    cd BD_BWT_index
    cmake .
    make
    cd ..
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
