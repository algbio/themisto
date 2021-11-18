#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <set>
#include <algorithm>
#include "throwing_streams.hh"
#include "sdsl/wavelet_trees.hpp"
#include "WheelerIndex.hh"

using namespace std;

// Nodes are k-mers and edges are (k+1)-mers.
template<typename bitvector_t = sdsl::bit_vector, // this should be one of the sdsl bit vector types, or compatible.
         typename rank_select_string_t = sdsl::wt_huff<bitvector_t>> // this should be one of the sdsl wavelet tree types, or compatible.
class BOSS : public wgi::WheelerIndex<bitvector_t, rank_select_string_t>{

private:

    typedef wgi::WheelerIndex<bitvector_t, rank_select_string_t> Base;

    int64_t k;

public:

    // The smallest k we allow is k = 1. In that case nodes are 1-mers and edges 2-mers.
    // Value k=0 would be weird because then edges are just self-loops from the empty string.
    // But then we can ever add only one edge because a second one would violate the Wheeler
    // graph rule that all incoming edges to a node must have the same label.
    // However, a 1-mer as a dummy edgemer is still a thing. It represents and edge from
    // the source node. 0-mer dummies don't make sense.
    BOSS() : k(1) {}
    BOSS(int64_t k) : k(k) {
        assert(k >= 1);
    }

    BOSS(const string& outlabels_string, const sdsl::bit_vector& indegs, const sdsl::bit_vector& outdegs, const vector<int64_t >& C, int64_t k) 
    : Base(outlabels_string, indegs, outdegs, C), k(k){
    }

    BOSS(const string& outlabels_string, const vector<bool>& indegs, const vector<bool>& outdegs, int64_t k) {
        sdsl::bit_vector sdsl_indegs(indegs.size()), sdsl_outdegs(outdegs.size());
        for(LL i = 0; i < (LL)indegs.size(); i++) sdsl_indegs[i] = indegs[i];
        for(LL i = 0; i < (LL)outdegs.size(); i++) sdsl_outdegs[i] = outdegs[i];

        vector<LL> counts(256);
        for(char c : outlabels_string) counts[c]++;
        vector<LL> C = char_counts_to_C_array(counts);
        BOSS<bitvector_t> newboss(outlabels_string, sdsl_indegs, sdsl_outdegs, C, k);
        *this = newboss;
    }

    int64_t get_k() const {return k;}

    LL serialize(ostream& os) const{
        LL written = 0;
        written += Base::serialize(os);
        os.write((char*)&k, sizeof(k));
        written += sizeof(k);
        return written;
    }

    void load(istream& is){
        Base::load(is);
        is.read((char*)&k, sizeof(k));
    }

    // Returns the node id of the given k-mer, or -1 if it does not exist in the BOSS.
    int64_t find_kmer(const char* kmer) const{
        // [l,r] = current colexicographic range (zero-indexed)
        int64_t l, r;
        std::tie(l,r) = Base::search(kmer, k);
        if(r < l) return -1;
        return l;        
    }

    int64_t find_kmer(const string& S) const{
        assert((LL)S.size() == k);
        return find_kmer(S.c_str());
    }

    // Takes as input a string S of length m with 1 <= m <= k+1.
    // Returns whether there is an edge from the node of S[0..m-2] 
    // to the node of S[1..m-1].
    bool edgemer_exists(const string& edgemer) const{
        return edgemer_exists(edgemer.c_str(), edgemer.size());
    }

    // Takes len such that 1 <= len <= k+1.
    // Returns whether there is an edge from the node of S[0..len-2] 
    // to the node of S[1..len-1].
    // If len < k+1, then the search is done from the root.
    bool edgemer_exists(const char* S, int64_t len) const{
        assert(len <= k+1);
        assert(len >= 1); // Edgemers have length at least 2, except dummmy edgemers can be 1-mers. But not 0-mers.
        return get_edgemer_destination(S,len) != -1;
    }

    // Takes len such that 1 <= len <= k+1.
    // Returns the node at the end of a path labeled with S.
    // If does not exist, returns -1.
    // If len < k+1, then the search is done from the root.
    int64_t get_edgemer_destination(const char* S, int64_t len) const{
        assert(len <= k+1);
        assert(len >= 1); // Edgemers have length at least 2, except dummmy edgemers can be 1-mers. But not 0-mers.
        if(len == k+1){
            int64_t node = find_kmer(S); // Find the node of the prefix
            return walk(node, S[len-1]); // Take the last step
        } else{
            return search_from_root(S, len);
        }
    }

    // Computes the label of the given node to the give char array. The char array label should
    // have at least k+1 bytes of space: k bytes for the k-mer and 1 byte for the nuint64_t terminator.
    // The label may be shorter than k characters. Returns the length of the label.
    int64_t get_node_label(int64_t node, char* label) const{
        assert(node != -1);

        label[k] = '\0';
        for(int64_t i = 0; i < k; i++){       
            // Basically do what we do in the forward walk function, but in reverse
            int64_t pos_in_indegs = Base::indegs_select1(node+1);
            if(pos_in_indegs == Base::indegs_size()-1 || Base::indegs_at(pos_in_indegs+1) == 1){
                // Indegree zero. Shift the current label to the start and return.
                int64_t old_start = k-1-(i-1);
                for(int64_t j = 0; j < i; j++)
                    label[j] = label[old_start + j];
                label[i] = '\0';
                return i;
            } 
            int64_t edge_rank = pos_in_indegs - node; // zero-based rank

            // There are edge_rank edges before this edge. We want the largest c such
            // that C[c] <= edge_rank. This c is the incoming label to the node
            char c = 0;
            for(int64_t char_idx = 0; char_idx < Base::alphabet_size(); char_idx++){
                char a = Base::alphabet_at(char_idx);
                if(Base::C_array_at(a) <= edge_rank) c = a;
            }
            label[k-1-i] = c;

            int64_t rank_within_c_block = edge_rank - Base::C_array_at(c); // zero-based rank
            int64_t outedge_pos_in_outlabels = Base::outlabels_select(rank_within_c_block+1, c);
            int64_t outedge_pos_in_outdegs = Base::outdegs_select0(outedge_pos_in_outlabels+1);
            node = Base::outdegs_rank1(outedge_pos_in_outdegs) - 1;
        }
        return k;
    }

    // todo: test
    string get_node_label(int64_t node) const{
        char label[k+1];
        get_node_label(node, label);
        string S(label);
        return S;
    }

    /**
     * \param node
     * \param c
     * \return A node that is reached by following an edge labeled with c from node. If such node does not exist, returns -1. If node is -1, returns -1.
     */
    int64_t walk(int64_t node, char c) const{

        if(node == -1) return -1;

        int64_t start, end;
        std::tie(start,end) = Base::outedge_range(node);
        if(end < start) return -1;

        // The labels of outgoing nodes are in outlabels[start..end]

        int64_t c_start_rank = Base::outlabels_rank(start, c);
        int64_t c_end_rank = Base::outlabels_rank(end+1, c);
        if(c_start_rank == c_end_rank) return -1; // Edge not found
        
        int64_t edge_rank = Base::C_array_at(c) + c_end_rank; // \in [1, number of edges]
        int64_t pos_in_indegs = Base::indegs_select0(edge_rank);
        int64_t destination_rank = Base::indegs_rank1(pos_in_indegs); // \in [1, number of nodes]
        return destination_rank - 1; // -1: back to zero indexing of nodes
    }


    /**
     * Returns the node reached by walking from the source node (indegree 0) of the
     * graph following the given path label.
     * 
     * \param pathlabel The search string
     * \param length The length of the search string
     * \return The id of the node at the end of the path, or -1 if does not exist.
     */
    int64_t search_from_root(const char* pathlabel, int64_t length) const{
        int64_t node = 0;
        for(int64_t i = 0; i < length; i++){
            node = walk(node, pathlabel[i]);
            if(node == -1) return -1;
        }
        return node;
    }

    /**
     * Returns the node reached by walking from the source node (indegree 0) of the
     * graph following the given path label.
     * 
     * \param pathlabel The search string
     * \return The id of the node at the end of the path, or -1 if does not exist.
     */    
    int64_t search_from_root(const string& pathlabel) const{
        return search_from_root(pathlabel.c_str(), pathlabel.size());
    }

    /**
     *  Returns all edgemers of length (k+1).
     */
    set<string> get_full_edgemers() const{
        set<string> edgemers;
        for(LL node = 0; node < Base::number_of_nodes(); node++){
            string label = get_node_label(node);
            if((LL)label.size() == k){
                string out = Base::node_outlabels(node);
                for(char c : out){
                    edgemers.insert(label + c);
                }
            }
        }
        return edgemers;
    }

    // Assumes edgemer exists
    LL edgemer_string_to_edge_wheeler_rank(const string& edgemer) const{
        LL node = find_kmer(edgemer.c_str()); // Searches first k chars
        assert(node != -1);
        char c = edgemer.back();
        LL edge_idx = Base::outdegs_rank0(Base::outdegs_select1(node+1));
        return Base::C_array_at(c) + Base::outlabels_rank(edge_idx, c);
    }

    bool has_exactly_one_source_node() const{
        if(Base::indegree(0) != 0) return false; // Node 0 must be a source node
        for(LL node = 1; node < Base::number_of_nodes(); node++){
            if(Base::indegree(node) == 0) return false; // Node >= 1 must not be a source node
        }
        return true;
    }

    vector<bool> get_dummy_node_marks() const{
        LL count = 0;
        vector<pair<LL,LL>> dfs_stack; // pairs (node, depth)
        dfs_stack.push_back({0, 0});
        // dfs to depth k-1
        // the dummy part is a tree so no visited-list is required

        vector<bool> marks(Base::number_of_nodes());
        LL v,d; // node,depth
        while(!dfs_stack.empty()){
            tie(v,d) = dfs_stack.back();
            dfs_stack.pop_back();
            if(d < k){
                count++;
                marks[v] = 1;
            }
            if(d < k-1){ // Push children
                LL out_l, out_r; tie(out_l, out_r) = Base::outedge_range(v);
                for(LL i = out_l; i <= out_r; i++){
                    LL edge_wr = Base::outedge_index_to_wheeler_rank(i);
                    dfs_stack.push_back({Base::edge_destination(edge_wr), d+1});
                }
            }
        }
        return marks;
    }

    vector<LL> char_counts_to_C_array(const vector<LL>& counts){
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

    // For debugging.
    set<string> get_all_edgemers() const{
        set<string> ans;
        for(LL v = 0; v < Base::number_of_nodes(); v++){
            string x = get_node_label(v);
            if(x.size() == get_k()){
                for(char c : Base::node_outlabels(v))
                    ans.insert(x + c);
            }
        }
        return ans;
    }

    template<typename T>
    friend bool operator==(BOSS<T>& boss1, BOSS<T>& boss2);

};

template<typename bitvector_t>
bool operator==(BOSS<bitvector_t>& boss1, BOSS<bitvector_t>& boss2){
    return *static_cast<wgi::WheelerIndex<bitvector_t>*>(&boss1) == 
           *static_cast<wgi::WheelerIndex<bitvector_t>*>(&boss2) && boss1.k == boss2.k;
}