#pragma once

#include <string>
#include "sbwt/SBWT.hh"
#include "new_coloring.hh"
#include "SeqIO.hh"
#include "variants.hh"

using namespace std;
using namespace sbwt;

typedef int64_t color_t;

void pseudoalign(const plain_matrix_sbwt_t& SBWT, const Coloring<>& coloring, int64_t n_threads, std::string inputfile, std::string outputfile, bool reverse_complements, int64_t buffer_size, bool gzipped, bool sorted_output);

// returns a vector where element i is the ref ids aligned with query i
vector<set<LL> > parse_pseudoalignment_output_format_from_disk(string filename);
