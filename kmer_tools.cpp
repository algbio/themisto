/*#include "kmer_tools.hh"
#include "input_reading.hh"
#include "globals.hh"
#include "test_tools.hh"

#include <string>
#include <set>

using namespace std;

// Input and output specs in the header
void list_all_distinct_kmers(string input_fastafile, string outputfile, int64_t k){
    string concat;
    FASTA_reader fr(input_fastafile);
    while(!fr.done()){
        string read = fr.get_next_query_stream().get_all();
        concat += read + read_separator;
    }
    concat += (char)0x01;
    set<string> kmers_set = get_all_distinct_cyclic_kmers(concat, k);
    ofstream out(outputfile);
    for(string kmer : kmers_set){
        out << kmer << "\n";
    }
    out.close();
}

// Input and output specs in the header
void colex_sort_kmers(string kmerfile, string outfile){
    string line;
    ifstream in(kmerfile);
    vector<string> kmers;
    while(getline(in,line)){
        kmers.push_back(line);
    }
    in.close();

    sort(kmers.begin(), kmers.end(), colex_compare);

    ofstream out(outfile);
    for(string kmer : kmers){
        out << kmer << "\n";
    }
    out.close();
}*/