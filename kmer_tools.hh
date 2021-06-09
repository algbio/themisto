#pragma once
#include <string>
#include "input_reading.hh"
#include "globals.hh"
#include "test_tools.hh"
#include "KMC_wrapper.hh"
#include <set>

using namespace std;

// Takes and output file from list_all_distinct_kmers
// Writes into outfile the lines of the input file sorted
// in colexicographic order
void colex_sort_kmers_in_memory(string kmerfile, string outfile){
    string line;
    throwing_ifstream in(kmerfile);
    vector<string> kmers;
    while(in.getline(line)){
        kmers.push_back(line);
    }
    in.close();

    sort(kmers.begin(), kmers.end(), colex_compare);

    throwing_ofstream out(outfile);
    for(string kmer : kmers){
        out << kmer << "\n";
    }
    out.close();
}
