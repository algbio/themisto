# About Themisto
Themisto is a succinct colored de Bruijn graph supporting pseudo-alignment against a database of reference sequences similar to the tool Kallisto. For more information, see the [webpage](https://www.helsinki.fi/en/researchgroups/genome-scale-algorithmics/themisto) and the [paper](https://www.biorxiv.org/content/biorxiv/early/2020/04/04/2020.04.03.021501/DC1/embed/media-1.pdf?download=true). The pseudoalignment algorithm is modeled after the tool Kallisto.

## Colored de Bruijn graph definition

The de Bruijn graph is defined so that the nodes represent k-mers and edges (k+1)-mers. The graph is [edge centric](https://www.biostars.org/p/175058/#256741), meaning that there is an edge from u to v if there is a (k+1)-mer in the data that is suffixed by u and prefixed by v. The set of nodes is the set of endpoints of all edges. Note that this implies that orphan k-mers (those that are not connected to any edge) are not in the graph. Each edge is associated with a set of colors. The color set contains the colors of all input sequences that have the (k+1)-mer corresponding to the edge. Each node is also given a color set that is the union of the color sets of the edges connected to the node.

We use the KMC3 library to list the distinct (k+1)-mers to construct the graph. Since KMC3 only works with the DNA alphabet ACGT, we must preprocess the data so that it does not have any characters outside of the alphabet ACGT. By default, we delete all (k+1)-mers that contain a letter that is outside of the alphabet ACGT. We also offer an option `--randomize-non-ACGT` that replaces non-ACGT characters with random nucleotides. If you would like to deal with non-ACGT characters differently, please preprocess the data before running.

# Installation
## Requirements

We currently support only Linux and macOS. For compilation, you will need a C++17 compliant compiler with OpenMP support, and CMake v3.1 or newer. If compiling with g++, make sure that the version is at least g++-8, or you might run into compilation errors with the standard library &lt;filesystem&gt; header.

## Compiling

### Linux

These instructions have been tested to work on Ubuntu 18.04. To build the software, enter the Themisto directory and run

```
cd build
cmake .. -DMAX_KMER_LENGTH=31 -DCMAKE_BUILD_ZLIB=1 -DCMAKE_BUILD_BZIP2=1
make
```

Where 31 is the maximum k-mer length (node length) to support, up to 255. The larger the k-mer length, the more time and memory the index construction takes. Values that are one less than a multiple of 32 work the best. This will create the binary at`build/bin/themisto`.

**Troubleshooting**: If you run into problems involving the &lt;filesystem&gt; header, you probably need to update your compiler. The compiler `g++-8` should be sufficient. Install a new compiler and direct CMake to use it with the `-DCMAKE_CXX_COMPILER` option. For example, to set the compiler to `g++-8`, run CMake with the option `-DCMAKE_CXX_COMPILER=g++-8`. 

## MacOS

Compiling Themisto on macOS requires users to first install gcc-8 or newer from homebrew with

```
brew install gcc@8
```

Afterwards, Themisto can be compiled by entering the directory and running

```
cd build
cmake -DCMAKE_C_COMPILER=$(which gcc-8) -DCMAKE_CXX_COMPILER=$(which g++-8) -DMAX_KMER_LENGTH=31 ..
make
```

Where 31 is the maximum k-mer length to support, up to 255. The larger the k-mer length, the more time and memory the index construction takes. If you run into problems involving zlib or bzip2, add `-DCMAKE_BUILD_ZLIB=1 -DCMAKE_BUILD_BZIP2=1` into the cmake command.

Note that macOS has a very small limit for the number of concurrently opened files. Themisto can use temporary files to conserve RAM, and may run into this limit. To increase the limit, run the command

```
ulimit -n 2048
```

# Usage

## Quick start

Themisto takes as an input a set of sequences in FASTA or FASTQ format, and a file specifying the color (a non-negative integer) of each sequence. The i-th line of the color file contains the color of the i-th sequence in the sequence file. For optimal compression, use color numbers in the range [0, n-1], where n is the number of distinct colors. If no color file is given, the index is built without colors. This way, the user can later try multiple colorings without recomputing the de Bruijn graph.

There is an example dataset with sequences at `example_input/coli3.fna` and colors at `example_input/colors.txt`. To build the index with order k = 30, such that the index files are written to `my_index.themisto.dbg` and `my_index.themisto.colors`, using the directory `temp` as temporary storage, using four threads and up to 1GB of memory.

```
./build/bin/themisto build --node-length 30 -i example_input/coli3.fna -c example_input/colors.txt --index-prefix my_index --temp-dir temp --mem-megas 1000 --n-threads 4
```

We recommend to use a fast SSD drive for the temporary directory. With a reasonable desktop workstation and an SSD drive, the program should take about one minute on this example input. Beware: for inputs that are in the range of tens of gigabytes, the index construction may need over a terabyte of temporary disk space.

To align the four sequences in `example_input/queries.fna` against the index we just built, writing output to out.txt run:

```
./build/bin/themisto pseudoalign --query-file example_input/queries.fna --index-prefix my_index --temp-dir temp --out-file out.txt --n-threads 4
```

This should produce the following output file:

```
0 43 748 
1 524 
2 855 
3 787 
```

There is one line for each query sequence. The lines may appear in a different order if parallelism was used. The first integer on a line is the 0-based rank of a query sequence in the query file, and the rest of the integers are the colors that are pseudoaligned with the query. For example, here the query with rank 2 (i.e. the 3rd sequence in the query file) pseudoaligns to color 855.

## Full instructions for index construction

This command builds an index consisting of compact de Bruijn graph using the BOSS data structure (implemented as a [Wheeler graph](https://www.sciencedirect.com/science/article/pii/S0304397517305285)) and color information. The input is a set of reference sequences in a single file in fasta or fastq format, and a colorfile, which is a plain text file containing the colors (integers) of the reference sequences in the same order as they appear in the reference sequence file, one line per sequence.

```
Usage:
  build [OPTION...]

  -k, --node-length arg         The k of the k-mers.
  -i, --input-file arg          The input sequences in FASTA or FASTQ 
                                format. The format is inferred from the 
                                file extension. Recognized file extensions 
                                for fasta are: .fasta, .fna, .ffn, .faa and 
                                .frn . Recognized extensions for fastq are: 
                                .fastq and .fq . If the file ends with .gz, 
                                it is uncompressed into a temporary 
                                directory and the temporary file is deleted 
                                after use.
  -c, --color-file arg          One color per sequence in the fasta file, 
                                one color per line. If not given, the 
                                sequences are given colors 0,1,2... in the 
                                order they appear in the input file. 
                                (default: "")
  -o, --index-prefix arg        The de Bruijn graph will be written to 
                                [prefix].themisto.dbg and the color 
                                structure to [prefix].themisto.colors. If 
                                not given, the filename of the sequence 
                                file used as the prefix.
      --temp-dir arg            Directory for temporary files. This 
                                directory should have fast I/O operations 
                                and should have as much space as possible.
  -m, --mem-megas arg           Number of megabytes allowed for external 
                                memory algorithms. Default: 1000 (default: 
                                1000)
  -t, --n-threads arg           Number of parallel exectuion threads. 
                                Default: 1 (default: 1)
      --randomize-non-ACGT      Replace non-ACGT letters with random 
                                nucleotides. If this option is not given, 
                                (k+1)-mers containing a non-ACGT character 
                                are deleted instead.
  -d, --colorset-pointer-tradeoff arg
                                This option controls a time-space tradeoff 
                                for storing and querying color sets. If 
                                given a value d, we store color set 
                                pointers only for every d nodes on every 
                                unitig. The higher the value of d, the 
                                smaller then index, but the slower the 
                                queries. The savings might be significant 
                                if the number of distinct color sets is 
                                small and the graph is large and has long 
                                unitigs. (default: 1)
      --no-colors               Build only the de Bruijn graph without 
                                colors.
      --load-dbg                If given, loads a precomputed de Bruijn 
                                graph from the index directory. If this is 
                                given, the parameter -k must not be given 
                                because the order k is defined by the 
                                precomputed de Bruijn graph.
  -h, --help                    Print usage
```

## Full instructions for `pseudoalign`

This command aligns query sequences against an index that has been built previously. The output is one line per input read. Each line consists of a space-separated list of integers. The first integer specifies the rank of the read in the input file, and the rest of the integers are the identifiers of the colors of the sequences that the read pseudoaligns with. If the program is ran with more than one thread, the output lines are not necessarily in the same order as the reads in the input file. This can be fixed with the option --sort-output.

The query can be given as one file, or as a file with a list of files. In the former case, we must specify one output file with the options --out-file, and in the latter case, we must give a file that lists one output filename per line using the option --out-file-list.

The query file(s) should be in fasta of fastq format. The format is inferred from the file extension. Recognized file extensions for fasta are: .fasta, .fna, .ffn, .faa and .frn . Recognized extensions for fastq are: .fastq and .fq

```
Usage:
  pseudoalign [OPTION...]

  -q, --query-file arg       Input file of the query sequences (default: 
                             "")
      --query-file-list arg  A list of query filenames, one line per 
                             filename (default: "")
  -o, --out-file arg         Output filename. (default: "")
      --out-file-list arg    A file containing a list of output filenames, 
                             one per line. (default: "")
  -i, --index-prefix arg     The index prefix that was given to the build 
                             command.
      --temp-dir arg         Directory for temporary files.
      --rc                   Whether to to consider the reverse complement 
                             k-mers in the pseudoalignemt.
  -t, --n-threads arg        Number of parallel exectuion threads. Default: 
                             1 (default: 1)
      --gzip-output          Compress the output files with gzip.
      --sort-output          Sort the lines of the out files by sequence 
                             rank in the input files.
  -h, --help                 Print usage
```

Examples:

Pseudoalign reads.fna against an index:
```
./build/bin/themisto pseudoalign --query-file reads.fna --index-prefix my_index --temp-dir temp --out-file out.txt
```

Pseudoalign a list of fasta files in input_list.txt into output filenames in output_list.txt:
```
./build/bin/themisto pseudoalign --query-file-list input_list.txt --index-prefix my_index --temp-dir temp --out-file-list output_list.txt
```

Pseudoalign reads.fna against an index using also reverse complements:
```
./build/bin/themisto pseudoalign --rc --query-file reads.fna --index-prefix my_index --temp-dir temp --outfile out.txt
```

## Extracting unitigs with `extract-unitigs`

This program dumps the unitigs and optionally their colors out of an existing Themisto index.

```
Usage:
  extract-unitigs [OPTION...]

  -i, --index-prefix arg  The index prefix that was given to the build 
                          command.
      --unitigs-out arg   Output filename for the unitigs (outputted in 
                          FASTA format). (default: "")
      --colors-out arg    Output filename for the unitig colors. If this 
                          option is not given, the colors are not computed. 
                          Note that giving this option affects the unitigs 
                          written to unitigs-out: if a unitig has nodes 
                          with different color sets, the unitig is split 
                          into maximal segments of nodes that have equal 
                          color sets. The file format of the color file is 
                          as follows: there is one line for each unitig. 
                          The lines contain space-separated strings. The 
                          first string on a line is the FASTA header of a 
                          unitig, and the following strings on the line are 
                          the integer color labels of the colors of that 
                          unitig. The unitigs appear in the same order as 
                          in the FASTA file. (default: "")
  -h, --help              Print usage
```

# For developers: building the tests

```
git submodule init
git submodule update
cd googletest
mkdir build
cd build
cmake ..
make
cd ../../build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DMAX_KMER_LENGTH=255
make
```

This builds the tests to build/bin/themisto_tests

# License

This software is licensed under GPLv2. See LICENSE.txt.
