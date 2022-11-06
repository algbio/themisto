#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include <cassert>
#include <set>
#include <unordered_map>
#include <map>
#include "globals.hh"
#include "sbwt/globals.hh"
#include "sbwt/throwing_streams.hh"
#include "sdsl/io.hpp" // stdlib printing
#include "sbwt/variants.hh"
#include "sbwt/SBWT.hh"
#include <cassert>

using namespace std;

template <typename T>
void print(T C){
    cout << "[";
    for(auto x : C) cout << x << ", ";
    cout << "]" << endl;
}

template <typename T, typename ostream_t>
void print(T C, ostream_t& out){
    out << "[";
    for(auto x : C) out << x << ", ";
    out << "]" << endl;
}

template<typename T>
void write_vector(vector<T>& v, string outfilename){
    sbwt::throwing_ofstream out(outfilename);
    for(T& t : v) out << t << "\n";
}

set<char> get_alphabet(string S);

set<string> get_all_distinct_kmers(string S, int64_t k);

set<string> get_all_distinct_cyclic_kmers(string S, int64_t k);

set<string> get_all_distinct_cyclic_kmers(vector<string>& A, int64_t k);

vector<string> get_all_kmers(string& S, int64_t k);

vector<string> all_binary_strings_up_to(int64_t k);

string get_random_dna_string(int64_t length, int64_t alphabet_size);

string get_random_string(int64_t length, int64_t alphabet_size);

vector<string> get_sorted_suffixes(string S);

void write_as_fasta(vector<string>& seqs, string fasta_filename);
void write_as_fastq(vector<string>& seqs, string fastq_filename);

vector<string> dump_node_labels(sbwt::plain_matrix_sbwt_t& SBWT);

template<typename T>
T to_disk_and_back(T& c){
    string f = sbwt::get_temp_file_manager().create_filename();
    sbwt::throwing_ofstream out(f, ios::binary);
    c.serialize(out.stream);
    out.close();

    T c_loaded;
    sbwt::throwing_ifstream in(f, ios::binary);
    c_loaded.load(in.stream);
    return c_loaded;
}
