#pragma once

#include <string>
#include "sbwt/SBWT.hh"
#include "new_coloring.hh"
#include "SeqIO.hh"
#include "variants.hh"

using namespace std;
using namespace sbwt;

typedef uint64_t color_t;

// Stores the intersection into buf1 and returns the number of elements in the
// intersection (does not resize buf1). Buffer elements must be sorted.
// Assumes all elements in a buffer are distinct
LL intersect_buffers(vector<color_t>& buf1, LL buf1_len, vector<color_t>& buf2, LL buf2_len);

// Stores the union into result_buf and returns the number of elements in the
// union (does not resize result_buf). Buffers elements must be sorted.
// Assumes all elements in a buffer are distinct. result_buf must have enough
// space to accommodate the union
LL union_buffers(vector<color_t>& buf1, LL buf1_len, vector<color_t>& buf2, LL buf2_len, vector<color_t>& result_buf);

void pseudoalign(const plain_matrix_sbwt_t& SBWT, const Coloring<>& coloring, int64_t n_threads, std::string inputfile, std::string outputfile, bool reverse_complements, int64_t buffer_size, bool gzipped, bool sorted_output);

// returns a vector where element i is the ref ids aligned with query i
vector<set<LL> > parse_pseudoalignment_output_format_from_disk(string filename);
