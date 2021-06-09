#pragma once

#include <string>

using namespace std;

// Only works for alphabet {a,c,g,t,A,C,G,T}
// Return the prefix of the databse filename
extern string KMC_wrapper(int64_t k, int64_t ram_gigas, int64_t n_threads, string fastafile, string tempdir, bool only_canonical_kmers);
