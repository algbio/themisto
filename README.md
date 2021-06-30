# About Themisto
Themisto is a compact colored de Bruijn graph supporting pseudo-alignment against a database of reference sequences similar to the tool Kallisto. For more information, see https://www.helsinki.fi/en/researchgroups/genome-scale-algorithmics/themisto.

# Installation
## Requirements
Compilation: C++17 compliant compiler with OpenMP support, and CMake v3.1 or newer.

## Compiling
Enter the Themisto directory and run


```
	cd build
	cmake .. -DMAX_KMER_LENGTH=60
    make
```

Where 60 is the maximum k-mer length to support, up to 255. The larger the k-mer length, the more time and memory the index construction takes.

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
Examples:

Build BOSS and colors:
```

  ./build/bin/build_index -k 31 --mem-megas 10000 --input-file references.fna --color-file colors.txt --index-dir index --temp-dir temp
```

Build only the BOSS
```
  ./build/bin/build_index -k 31 --mem-megas 10000 --input-file references.fna --index-dir index --temp-dir temp
```

Load a previously built BOSS from the index directory and compute the colors:
```
  ./build/bin/build_index --mem-megas 10000 --input-file references.fna --color-file colors.txt --index-dir index --temp-dir temp --load-boss
```

The program will then print instructions on how to use it.

## Pseudoalignment
To pseudoalign against the index, run:

    pseudoalign

The program will then print instructions on how to use it.

# License

This software is licensed under GPLv2. See LICENSE.txt.
