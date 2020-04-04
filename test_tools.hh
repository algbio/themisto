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
#include <cassert>

using namespace std;
typedef int64_t LL;

template <typename T>
void print(T C){
    cout << "[";
    for(auto x : C) cout << x << ", ";
    cout << "]" << endl;
}

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

template <typename S, typename T>
ostream& operator<<(ostream& os, const unordered_map<S,T>& v){
    os << "[";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << it->first << ": " << it->second;
    }
    os << "]";
    return os;
}

template <typename S, typename T>
ostream& operator<<(ostream& os, const map<S,T>& v){
    os << "{";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << it->first << ": " << it->second;
    }
    os << "}";
    return os;
}

template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v){
    os << "[";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << *it;
    }
    os << "]";
    return os;
}

template <typename T>
ostream& operator<<(ostream& os, const set<T>& v){
    os << "[";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << *it;
    }
    os << "]";
    return os;
}

template <typename T>
ostream& operator<<(ostream& os, const multiset<T>& v){
    os << "[";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << *it;
    }
    os << "]";
    return os;
}

template <typename S, typename T>
ostream& operator<<(ostream& os, const pair<S,T>& x){
    os << "(" << x.first << ", " << x.second << ")";
    return os;
}