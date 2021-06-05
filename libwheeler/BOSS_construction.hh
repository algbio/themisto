#pragma once
#include "BOSS.hh"
//#include "util/globals.hh"
//#include "Kmer.hh"

vector<LL> char_counts_to_C_array(vector<LL> counts){
    vector<LL> C(256); // Cumulative sum of counts

    // Compute cumulative sum of counts
    for(LL i = 0; i < (LL)C.size(); i++){
        C[i] = counts[i];
        if(i > 0) C[i] += C[i-1];
    }

    // Shift C to the right by one because that's how it's defined
    for(LL i = 256-1; i >= 0; i--){
        if(i == 0) C[i] = 0;
        else C[i] = C[i-1];
    }

    return C;

}

// Constructs a BOSS where nodes are k-mers and two nodes are connected by an
// edge iff the (k+1) that is formed by merging the k-mers exists in the reads
template<typename bitvector_t = sdsl::bit_vector>
BOSS<bitvector_t> build_BOSS_with_maps(const vector<string>& reads, LL k){
    map<string, LL> k_counts; // k-mer counts. Keys are reverses of the k-mers so that the order is colexicographic
    map<string, LL> k_plus_1_counts; // k+1 mer counts. Keys are reverses of the k-mers so that the order is colexicographic

    set<char> alphabet;
    for(string S : reads) for(char c : S) alphabet.insert(c);
    alphabet.insert('$');

    string dummy_prefix(k,'$');
    for(string S : reads){
        S = dummy_prefix + S;
        std::reverse(S.begin(), S.end());
        for(LL i = 0; i < (LL)S.size() - k + 1; i++)
            k_counts[S.substr(i,k)]++;
        for(LL i = 0; i < (LL)S.size() - k; i++)
            k_plus_1_counts[S.substr(i,k+1)]++;
    }

    string outlabels;
    vector<bool> outdegrees; // Concatenation of unary representations. Indegree k -> 1 0^k
    vector<bool> indegrees;

    for(auto P : k_counts){
        // In colex order
        string rev_kmer = P.first;
        LL indegree = 0;
        LL outdegree = 0;
        for(char c : alphabet){
            if(k_plus_1_counts[c + rev_kmer] > 0){
                outlabels += c;
                outdegree++;
            }
            if(k_plus_1_counts[rev_kmer + c] > 0){
                indegree++;
            }
        }

        outdegrees.push_back(1);
        for(LL i = 0; i < outdegree; i++){
            outdegrees.push_back(0);
        }

        indegrees.push_back(1);
        for(LL i = 0; i < indegree; i++){
            indegrees.push_back(0);
        }
    }

    sdsl::bit_vector sdsl_indegs(indegrees.size()), sdsl_outdegs(outdegrees.size());
    for(LL i = 0; i < (LL)indegrees.size(); i++) sdsl_indegs[i] = indegrees[i];
    for(LL i = 0; i < (LL)outdegrees.size(); i++) sdsl_outdegs[i] = outdegrees[i];

    vector<LL> counts(256);
    for(char c : outlabels) counts[c]++;
    vector<LL> C = char_counts_to_C_array(counts);
    return BOSS<bitvector_t>(outlabels, sdsl_indegs, sdsl_outdegs, C, k);
}