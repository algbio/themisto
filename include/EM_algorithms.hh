#pragma once

#include <string>
#include <vector>
#include <fstream>
#include "globals.hh"
#include "EM_sort.hh"

// Returns output filename
string EM_sort_big_endian_LL_pairs(string infile, LL ram_bytes, LL key, LL n_threads);

// Returns output filename
string EM_delete_duplicate_LL_pair_records(string infile);

// (node, color) pairs to (node, color list) pairs
// Returns output filename
string EM_collect_colorsets_binary(string infile);

// Input: records (node, color set), in the format where
// first we have the length of the whole record in bytes
// Returns output filename
string EM_sort_by_colorsets_binary(string infile, LL ram_bytes, LL n_threads);

// Returns output filename
string EM_collect_nodes_by_colorset_binary(string infile);
