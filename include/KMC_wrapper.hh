#pragma once

#include <string>

using namespace std;

// Only works for alphabet {a,c,g,t,A,C,G,T}
//void KMC_wrapper(int64_t k, int64_t ram_gigas, int64_t n_threads, string fastafile, string tempdir, string database_filename_prefix, bool only_canonical_kmers, bool silent);


// Returns the KMC database prefix and the number of distinct k-mers that had abundance within the given bounds
pair<string, int64_t> run_kmc(const vector<string>& input_files, LL k, LL n_threads, LL ram_gigas, int64_t min_abundance, int64_t max_abundance);