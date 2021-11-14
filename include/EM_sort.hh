#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <set>
#include <utility>
#include <tuple>
#include <cstring>
#include <cstdio>
#include <cassert>
#include "globals.hh"
#include "Block.hh"
#include "ParallelBoundedQueue.hh"
#include "generic_EM_classes.hh"

/*
Design:
  
There are two operating modes: EM_LINES and EM_BINARY. The mode determines
how records are represented in disk and memory.

EM_LINES:
  - The input is a file with one line per record. user-provided comparison function 
    takes as a parameter two NULL-terminated C-style strings which are lines in the
    input file.
  - Internal workings: The records are read from disk into an object of class Block. The 
    data of the block is stored in one large char-array. All lines are concatenated, 
    but newlines are replaced with NULL-terminator so that strcmp can be used easily
    for doing the comparisons if wanted. A vector<LL> 'starts' stores the 
    start positions of each record. Sorting is done by permuting 'starts' and keeping the 
    data in place. When the data is written back to diks, the null-terminators are again 
    replaced by newline-characters, so the resulting file is a permutation of the lines 
    of the original file.

EM_BINARY:
  - Records are in binary such that the first 8 bytes are the big-endian representation
    of a 64-bit integer x. Then follows x-8 bytes of data. The comparison function now takes
    as a parameter the pointers to the starts of the x byte block of data of the records.
    When passing around records, the size x is always at the start of the data array.
  - Internal workings: The records are read from disk into an object of class Block, which
    has a char-array which will store the concatenation of all records. A vector<LL> 'starts'
    stores the start position of each record. Sorting is done by permuting 'starts' and keeping
    the data in place. The data is written back to disk according to the permutation of starts.
*/

using namespace std;

static const int EM_LINES = 0;
static const int EM_VARIABLE_BINARY = 1;
static const int EM_CONSTANT_BINARY = 2;

// Interprets the strings as integers (no leading zeros allowed) and returns:
//     -1 if x < y
//      0 if x = y
//      1 if x > y
int compare_as_numbers(const char* x, const char* y);
bool memcmp_variable_binary_records(const char* x, const char* y);
void copy_file(string infile, string outfile, LL buf_size);

void merge_files_generic(const std::function<bool(const char* x, const char* y)>& cmp, LL& merge_count, Generic_Record_Reader* reader, Generic_Record_Writer* writer);

void EM_sort_generic(string infile, string outfile, const std::function<bool(const char* x, const char* y)>& cmp, LL RAM_bytes, LL k, Generic_Block_Producer* producer, vector<Generic_Block_Consumer*> consumers, Generic_Record_Reader* reader, Generic_Record_Writer* writer);

// Constant size records of record_size bytes each
void EM_sort_constant_binary(string infile, string outfile, const std::function<bool(const char* x, const char* y)>& cmp, LL RAM_bytes, LL k, LL record_size, LL n_threads);

// Binary format of record: first 8 bytes give the length of the record, then comes the record
// Line format of record: A char sequence terminated by \n\0. For the cmp function the terminator is just \0.
// k = k-way merge parameter
// mode = EM_LINES or mode = EM_VARIABLE_BINARY or mode = EM_CONSTANT_BINARY
void EM_sort(string infile, string outfile, const std::function<bool(const char* x, const char* y)>& cmp, LL RAM_bytes, LL k, LL n_threads, LL mode);