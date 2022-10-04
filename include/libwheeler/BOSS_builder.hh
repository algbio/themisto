#pragma once

#include "Kmer.hh"
#include "BOSS.hh"
#include "test_tools.hh"
#include <unordered_map>
#include <stdexcept>


// Constructs the edge-centric de-bruijn graph where nodes are k-mers.
// i.e. for every (k+1)-mer, adds and edge between the prefix and suffix k-mer.
// The set of nodes in implicitly defined as the set of endpoints of all these
// nodes.
// Warning: this is very inefficiently implemented! Used for debug purposes
template<typename bitvector_t = sdsl::bit_vector>
BOSS<bitvector_t> build_BOSS_with_maps(vector<string> reads, LL k, bool include_reverse_complements){

    if(include_reverse_complements){
        LL n = reads.size();
        for(LL i = 0; i < n; i++){
            reads.push_back(get_rc(reads[i]));
        }
    }

    struct Edge_Info{
        set<char> inlabels;
        set<char> outlabels;
    };

    struct colex_compare_object{
        bool operator()(const string& A, const string& B) const{
            return colex_compare(A,B);
        }
    };

    map<string, Edge_Info, colex_compare_object> M; // edgemer -> edge info

    // Add all edges
    for(string seq : reads){
        if(seq.size() >= k+1){
            for(string x : get_all_distinct_kmers(seq, k+1)){
                M[x.substr(0,k)].outlabels.insert(x.back());
                M[x.substr(1,k)].inlabels.insert(x[0]);
            }
        }
    }


    // Add dummy nodes
    map<string, Edge_Info, colex_compare_object> M_copy = M; // Take a copy so we don't edit M while we iterate it
    for(auto P : M_copy){
        string kmer = P.first;
        if(P.second.inlabels.size() == 0){
            // Need to add the prefixes
            for(LL len = 0; len <= k; len++){
                string prefix = kmer.substr(0,len);
                if(len != 0) M[prefix].inlabels.insert('$');
                if(len != k) M[prefix].outlabels.insert(kmer[len]);
            }
        }
    }

    // Make sure the root node exists
    if(M.find("") == M.end()) M[""] = Edge_Info();

    write_log("map build added " + to_string(M.size() - M_copy.size()) + " dummies (k = " + to_string(k) + ")", LogLevel::MAJOR);


    // Build the boss data structures    
    string outlabels; // Concatenation in outedge label sets in colexicographic order
    vector<bool> outdegrees; // Concatenation of unary representations. Indegree d -> 1 0^d
    vector<bool> indegrees; // Concatenation of unary representations. Outdegree d -> 1 0^d

    for(auto P : M){
        outdegrees.push_back(1);
        string debug_outlabels;
        for(char c : P.second.outlabels){
            outlabels.push_back(c);
            outdegrees.push_back(0);
        }

        indegrees.push_back(1);
        for(char c : P.second.inlabels){
            (void)c; // Unsed variable. Make compiler happy
            if(P.first != "") indegrees.push_back(0);
        }         
    }

    return BOSS<bitvector_t>(outlabels, indegrees, outdegrees, k);
}