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
#include "test_tools.hh"
#include <cassert>

using namespace std;
typedef int64_t LL;

set<char> get_alphabet(string S){
    set<char> ans;
    for(char c: S) ans.insert(c);
    return ans;
}

set<string> get_all_distinct_kmers(string S, LL k){
    set<string> kmers;
    for(LL i = 0; i < S.size()-k+1; i++){
        kmers.insert(S.substr(i,k));
    }
    return kmers;
}

set<string> get_all_distinct_cyclic_kmers(string S, LL k){
    set<string> kmers;
    for(LL i = 0; i < S.size(); i++){
        string kmer;
        for(LL j = 0; j < k; j++){
            kmer += S[(i+j) % S.size()];
        }
        kmers.insert(kmer);
    }
    return kmers;
}


set<string> get_all_distinct_cyclic_kmers(vector<string>& A, LL k){
    string concat;
    for(string read : A){
        concat += read_separator + read;
    }
    concat += '\x01'; // bibwt end sentinel

    return get_all_distinct_cyclic_kmers(concat,k);
}

vector<string> get_all_kmers(string& S, LL k){
    vector<string> kmers;
    for(LL i = 0; i < S.size()-k+1; i++){
        kmers.push_back(S.substr(i,k));
    }
    return kmers;
}

vector<string> all_binary_strings_up_to(int64_t k){ // For testing
    vector<string> ans;
    for(int64_t length = 1; length <= k; length++){
        for(int64_t mask = 0; mask < (1 << length); mask++){
            string s = "";
            for(int64_t i = 0; i < length; i++){
                if(mask & (1 << i)) s += 'A';
                else s += 'C';
            }
            ans.push_back(s);
        }
    }
    return ans;
}

string get_random_dna_string(int64_t length, int64_t alphabet_size){ // For testing
    string s;
    assert(alphabet_size >= 1 && alphabet_size <= 4);
    char alphabet[4] = {'A','T','G','C'};
    for(int64_t i = 0; i < length; i++){
        s.push_back(alphabet[rand() % alphabet_size]);
    }
    return s;
}

string get_random_string(int64_t length, int64_t alphabet_size){ // For testing
    string s;
    for(int64_t i = 0; i < length; i++){
        LL r = rand() % alphabet_size;
        s += 'a' + r;
    }
    return s;
}

vector<string> get_sorted_suffixes(string S){
    vector<string> suffixes;
    for(int64_t i = 0; i < S.size(); i++){
        suffixes.push_back(S.substr(i));
    }
    sort(suffixes.begin(), suffixes.end());
    return suffixes;
}

void write_as_fasta(vector<string>& seqs, string fasta_filename){
    throwing_ofstream out(fasta_filename);
    for(string& S : seqs) out << ">\n" << S << "\n";
}