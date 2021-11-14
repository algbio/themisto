#pragma once

#include <string>
#include <vector>
#include <fstream>
#include "globals.hh"
#include "EM_sort.hh"

string EM_sort_big_endian_LL_pairs(string infile, LL ram_bytes, LL key, LL n_threads);

// Input: file with one ';'-separated pair of strings (x,y) on each line
//        cmp_1 compares x values and cmp_y compares y values. 
//        The comparison functions should return < 0, 0 or > 0 like strcmp.
// Output: a file with the lines sorted by x if key = 0, or by y if key == 1
//         if the key values are equal, the other elements are compared
string EM_line_sort_pairs(string infile, int (*cmp_1) (const char* x_1, const char* x_2), 
                                         int (*cmp_2) (const char* x_1, const char* x_2), 
                                         LL key, LL ram_bytes, LL n_threads);

string EM_delete_duplicate_lines_from_sorted_file(string infile);

string debug_binary_to_string_form(string infile);

string EM_delete_duplicate_LL_pair_records(string infile);

// (node, color) pairs to (node, color list) pairs
string EM_collect_colorsets_binary(string infile);

// Input: records (node, color set), in the format where
// first we have the length of the whole record in bytes
string EM_sort_by_colorsets_binary(string infile, LL ram_bytes, LL n_threads);

string EM_collect_nodes_by_colorset_binary(string infile);


// Input: a file with ';'-separated distinct string pairs (x,y),
//        sorted by x if key = 0, and by y if key = 1
// Output: 
// if key_pos == 0:
//   then pairs (x,Y), where Y is a space-separated list of y which appear with x
// if key_pos == 1:
//   then pairs (X,y), where X is a space-separated list of x which appear with y
string EM_collect(string infile, LL key_pos);