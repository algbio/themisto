#pragma once

#include <string>

using namespace std;

// Only works for alphabet {a,c,g,t,A,C,G,T}
void KMC_wrapper(int64_t k, int64_t ram_gigas, int64_t n_threads, string fastafile, string tempdir, string database_filename_prefix, bool only_canonical_kmers);
