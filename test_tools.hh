#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
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

set<string> get_concat_kmers(vector<string>& A , vector<string>& B, LL k, char separator){
    string concat;
    for(string read : A){
        concat += separator + read;
    }
    for(string read : B){
        concat += separator + read;
    }
    concat += '\x01'; // bibwt end sentinel

    set<string> kmers_concat;
    string cyclic_concat = concat + concat; // Get kmers from this so that we get the wrap around at the end like in the bwt
    for(LL i = 0; i < cyclic_concat.size()-k+1; i++) kmers_concat.insert(cyclic_concat.substr(i,k));
    return kmers_concat;
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
                if(mask & (1 << i)) s += 'a';
                else s += 'b';
            }
            ans.push_back(s);
        }
    }
    return ans;
}

string get_random_string(int64_t length, int64_t alphabet_size){ // For testing
    string s;
    for(int64_t i = 0; i < length; i++){
        s.push_back('a' + rand() % alphabet_size);
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

class Argv{ // Class for turning a vector<string> into char**
private:

    // Forbid copying the class because it wont work right
    Argv(Argv const& other);
    Argv& operator=(Argv const& other);

public:

    char** array = NULL;
    int64_t size = 0;

    Argv(vector<string> v){
        array = (char**)malloc(sizeof(char*) * v.size());
        // Copy contents of v into array
        for(int64_t i = 0; i < v.size(); i++){
            char* s = (char*)malloc(sizeof(char) * (v[i].size() + 1)); // +1: space for '\0' at the end
            for(int64_t j = 0; j < v[i].size(); j++){
                s[j] = v[i][j]; // Can't use strcpy because s.c_str() is const
            }
            s[v[i].size()] = '\0';
            array[i] = s;
        }
        size = v.size();
    }

    ~Argv(){
        for(int64_t i = 0; i < size; i++) free(array[i]);
        free(array);
    }

};

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