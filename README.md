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

### Examples

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

### Full instructions

This program builds an index consisting of compact de Bruijn graph using the BOSS data structure and color information. The input is a set of reference sequences in a single file in fasta or fastq format, and a colorfile, which is a plain text file containing the colors of the reference sequences in the same order as they appear in the reference sequence file, one line per sequence. The names are given as ASCII strings, but they should not contain whitespace characters.

```
Usage:
  ./build/bin/build_index [OPTION...]

      --load-boss               If given, loads a precomputed BOSS from the
                                index directory
  -k, arg                       The k of the k-mers. Required only if
                                --load-boss is not given
  -i, --input-file arg          The input sequences in FASTA or FASTQ format.
                                The format is inferred from the file
                                extension. Recognized file extensions for fasta are:
                                .fasta, .fna, .ffn, .faa and .frn . Recognized
                                extensions for fastq are: .fastq and .fq . If
                                the file ends with .gz, it is uncompressed
                                into a temporary directory and the temporary file
                                is deleted after use.
  -c, --color-file arg          One color per sequence in the fasta file, one
                                color name per line. Required only if you
                                want to build the colors. (default: )
      --auto-colors             Instead of a color file let the program
                                automatically give colors integer names (0,1,2...)
  -o, --index-dir arg           Directory where the index will be built.
                                Always required, directory must exist before
                                running.
  -d, --colorset-pointer-tradeoff arg
                                This option controls a time-space tradeoff
                                for storing and querying color sets. If given a
                                value d, we store color set pointers only for
                                every d nodes on every unitig. The higher the
                                value of d, the smaller then index, but the
                                slower the queries. The savings might be
                                significant if the number of distinct color sets is
                                small and the graph is large and has long
                                unitigs. (default: 1)
      --temp-dir arg            Temporary directory. Always required,
                                directory must exist before running.
  -m, --mem-megas arg           Number of megabytes allowed for external
                                memory algorithms. Default: 1000 (default: 1000)
  -t, --n-threads arg           Number of parallel exectuion threads.
                                Default: 1 (default: 1)
      --pp-buf-siz arg          Size of preprocessing buffer (in bytes) for
                                fixing alphabet (default: 4096)
  -h, --help                    Print usage
```

## Pseudoalignment

### Examples

Pseudoalign reads.fna against an index:
```
  ./pseudoalign --query-file reads.fna --index-dir index --temp-dir temp --out-file out.txt
```

Pseudoalign a list of fasta files in input_list.txt into output filenames in output_list.txt:
```
  ./pseudoalign --query-file-list input_list.txt --index-dir index --temp-dir temp --out-file-list output_list.txt
```

Pseudoalign reads.fna against an index using also reverse complements:
```
  ./pseudoalign --rc --query-file reads.fna --index-dir index --temp-dir temp --outfile out.txt
```

### Full instructions

This program aligns query sequences against an index that has been built previously. The output is one line per input read. Each line consists of a space-separated list of integers. The first integer specifies the rank of the read in the input file, and the rest of the integers are the identifiers of the colors of the sequences that the read pseudoaligns with. If the program is ran with more than one thread, the output lines are not necessarily in the same order as the reads in the input file. This can be fixed with the option --sort-output.

If the coloring data structure was built with the --color-file option, then the integer identifiers of the colors can be mapped back to the provided color names by parsing the file coloring-mapping-id_to_name in the index directory. This file contains as many lines as there are distinct colors, and each line contains two space-separated strings: the first is the integer identifier of a color, and the second is the corresponding color name. In case the --auto-colors option was used, the integer identifiers are always numbers [0..n-1], where n is the total number of reference sequences, and the identifiers are assigned in the same order as the reference sequences were given to build_index.

 The query can be given as one file, or as a file with a list of files. In the former case, we must specify one output file with the options --out-file, and in the latter case, we must give a file that lists one output filename per line using the option --out-file-list.

The query file(s) should be in fasta of fastq format. The format is inferred from the file extension. Recognized file extensions for fasta are: .fasta, .fna, .ffn, .faa and .frn . Recognized extensions for fastq are: .fastq and .fq

```
Usage:
  ./build/bin/pseudoalign [OPTION...]

  -q, --query-file arg       Input file of the query sequences (default: )
      --query-file-list arg  A list of query filenames, one line per filename
                             (default: )
  -o, --out-file arg         Output filename. (default: )
      --out-file-list arg    A file containing a list of output filenames,
                             one per line. (default: )
  -i, --index-dir arg        Directory where the index will be built. Always
                             required, directory must exist before running.
      --temp-dir arg         Temporary directory. Always required, directory
                             must exist before running.
      --rc                   Whether to to consider the reverse complement
                             k-mers in the pseudoalignemt.
  -t, --n-threads arg        Number of parallel exectuion threads. Default: 1
                             (default: 1)
      --gzip-output          Compress the output files with gzip.
      --sort-output          Sort the lines of the out files by sequence rank
                             in the input files.
  -h, --help                 Print usage
```

# License

This software is licensed under GPLv2. See LICENSE.txt.
