#pragma once

#include <string>

using namespace std;

extern void KMC_wrapper(int64_t k, int64_t ram_gigas, int64_t n_threads, string fastafile, string outfile, string tempdir);
