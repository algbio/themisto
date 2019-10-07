#pragma once
#include <string>
#include "input_reading.hh"
#include "globals.hh"
#include "test_tools.hh"
#include "KMC_wrapper.hh"
#include <set>

using namespace std;


// Adds cyclic kmers that cross some separator.
// APPENDS to outfile, does not overwrite
void add_extra_kmers(string fastafile, string outfile, LL k){
    ofstream kmerfile(outfile, std::ofstream::out | std::ofstream::app);
    FASTA_reader fr(fastafile);
    string concat;
    while(!fr.done()){
        string read = fr.get_next_query_stream().get_all();
        if(read.size() < 2*k){
            concat += read_separator + read;
        } else{
            concat += read_separator + read.substr(0,k) + read.substr(read.size()-k,k);
        }
    }
    concat += BD_BWT_index<>::END;
    for(string kmer : get_all_distinct_cyclic_kmers(concat,k)){
        bool good = false;
        for(char c : kmer) if(c == read_separator || c == BD_BWT_index<>::END){
            good = true;
        }
        if(good) kmerfile << kmer << "\n";
    }
}

// Let S_1, S_2, ..., S_n be the sequences in the input fasta file.
// Let C = $ S_1 $ S_2 $ ... S_n X 
// where $ is the sequence separator defined in globals.hh and X is the byte 0x01
// Writes into outputfile all distinct cyclic k-mers of C, one k-mer per line, in any order
void list_all_distinct_cyclic_kmers_in_memory(string input_fastafile, string outputfile, int64_t k){
    string concat;
    FASTA_reader fr(input_fastafile);
    while(!fr.done()){
        string read = fr.get_next_query_stream().get_all();
        concat += read_separator + read;
    }
    concat += (char)0x01;
    set<string> kmers_set = get_all_distinct_cyclic_kmers(concat, k);
    ofstream out(outputfile);
    for(string kmer : kmers_set){
        out << kmer << "\n";
    }
    out.close();
}

// Let S_1, S_2, ..., S_n be the sequences in the input fasta file.
// Let C = $ S_1 $ S_2 $ ... S_n X 
// where $ is the sequence separator defined in globals.hh and X is the byte 0x01
// Writes into outputfile all distinct cyclic k-mers of C, one k-mer per line, in any order
void list_all_distinct_cyclic_kmers_in_external_memory(string input_fastafile, string outputfile, int64_t k, int64_t ram_bytes, int64_t n_threads){
    KMC_wrapper(k, max(ram_bytes / (LL)1e9, (LL)1), n_threads, input_fastafile, outputfile, temp_file_manager.get_dir());
    add_extra_kmers(input_fastafile, outputfile, k);
}

// Takes and output file from list_all_distinct_kmers
// Writes into outfile the lines of the input file sorted
// in colexicographic order
void colex_sort_kmers_in_memory(string kmerfile, string outfile){
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
}
