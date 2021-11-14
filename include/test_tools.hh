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
#include "stdlib_printing.hh"
#include <cassert>

using namespace std;
typedef int64_t LL;

template <typename T>
void print(T C){
    cout << "[";
    for(auto x : C) cout << x << ", ";
    cout << "]" << endl;
}

set<char> get_alphabet(string S);

set<string> get_all_distinct_kmers(string S, LL k);

set<string> get_all_distinct_cyclic_kmers(string S, LL k);

set<string> get_all_distinct_cyclic_kmers(vector<string>& A, LL k);

vector<string> get_all_kmers(string& S, LL k);

vector<string> all_binary_strings_up_to(int64_t k);

string get_random_dna_string(int64_t length, int64_t alphabet_size);

string get_random_string(int64_t length, int64_t alphabet_size);

vector<string> get_sorted_suffixes(string S);

void write_as_fasta(vector<string>& seqs, string fasta_filename);