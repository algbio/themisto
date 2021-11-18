#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <set>
#include <algorithm>
#include "sdsl/wavelet_trees.hpp"
#include "IO.hh"
#include <type_traits>
#include "globals.hh"

namespace wgi{

using namespace std;

// Based on the Wheeler graph representation of the De Bruijn graph, see:
// Gagie, T., Manzini, G., & Sirén, J. (2017). Wheeler graphs: A framework for BWT-based data structures. Theoretical computer science, 698, 67-78.
template<typename bitvector_t = sdsl::bit_vector, // this should be one of the sdsl bit vector types, or compatible.
         typename rank_select_string_t = sdsl::wt_huff<bitvector_t>> // this should be one of the sdsl wavelet tree types, or compatible.
class WheelerIndex{

private:

    rank_select_string_t outlabels;
    bitvector_t indegs;
    bitvector_t outdegs;
    vector<int64_t> C; // Array of length 256. C[c] is the number of edges with label strictly less than c in the graph
    
    typename bitvector_t::rank_1_type indegs_rs; // rank support
    typename bitvector_t::select_1_type indegs_ss1; // select support
    typename bitvector_t::select_0_type indegs_ss0; // select support
    
    typename bitvector_t::rank_1_type outdegs_rs; // rank support
    typename bitvector_t::select_1_type outdegs_ss1; // select support
    typename bitvector_t::select_0_type outdegs_ss0; // select support

    vector<char> alphabet;

    int64_t n_nodes;
    
    void set_supports(){
        this->indegs_rs.set_vector(&(this->indegs)); 
        this->indegs_ss1.set_vector(&(this->indegs));
        this->indegs_ss0.set_vector(&(this->indegs));
        this->outdegs_rs.set_vector(&(this->outdegs)); 
        this->outdegs_ss1.set_vector(&(this->outdegs));
        this->outdegs_ss0.set_vector(&(this->outdegs));
    }

    // An empty Index contains one node that represent the empty string.
    inline static const string EMPTY_WI_OUTLABELS = "";
    inline static const sdsl::bit_vector EMPTY_WI_INDEGS = {1};
    inline static const sdsl::bit_vector EMPTY_WI_OUTDEGS = {1};
    inline static const vector<int64_t > EMPTY_WI_C = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // 256 zeros

public:

    /**
     * \return Number of nodes in the graph.
     */
    int64_t number_of_nodes() const {return n_nodes;}

    /**
     * \return Number of edges in the graph.
     */
    int64_t number_of_edges() const {return outlabels.size();}

    /**
     * \return The pair (0, number of nodes -1).
     */
    pair<int64_t ,int64_t > full_node_range() const {return {0,n_nodes-1};}

    /**
     * \return A copy of the indegrees bitvector.
     */
    bitvector_t get_indegs() const {return indegs;}

    /**
     * \return The size of the indegrees bitvector.
     */
    int64_t indegs_size() const {return indegs.size();}

    /**
     * \param i Index in the indegrees bitvector.
     * \return The bit at index i.
     */
    int64_t indegs_at(int64_t i) const {return indegs[i];}

    /**
     * \param i Index in the indegrees bitvector.
     * \return The number of ones in the indegrees vector in the index range [0,i).
     */
    int64_t indegs_rank1(int64_t i) const {return indegs_rs.rank(i);}

    /**
     * \param i Index in the indegrees bitvector.
     * \return The number of zeros in the indegrees vector in the index range [0,i).
     */
    int64_t indegs_rank0(int64_t i) const {return i - indegs_rs.rank(i);}

    /**
     * \param 1
     * \return The position of the i-th one in the indegrees vector.
     */
    int64_t indegs_select1(int64_t i) const {return indegs_ss1.select(i);}

    /**
     * \param 1
     * \return The position of the i-th zero in the indegrees vector.
     */    
    int64_t indegs_select0(int64_t i) const {return indegs_ss0.select(i);}

    /**
     * \return A copy of the outdegrees bitvector.
     */
    bitvector_t get_outdegs() const {return outdegs;}

    /**
     * \return The size of the outdegrees bitvector.
     */
    int64_t outdegs_size() const {return outdegs.size();}

    /**
     * \param i Index in the outdegrees bitvector.
     * \return The bit at index i.
     */
    int64_t outdegs_at(int64_t i) const {return outdegs[i];}

    /**
     * \param i Index in the outdegrees bitvector.
     * \return The number of ones in the outdegrees vector in the index range [0,i).
     */
    int64_t outdegs_rank1(int64_t i) const {return outdegs_rs.rank(i);}

    /**
     * \param i Index in the outdegrees bitvector.
     * \return The number of zeros in the outdegrees vector in the index range [0,i).
     */
    int64_t outdegs_rank0(int64_t i) const {return i - outdegs_rs.rank(i);}

    /**
     * \param i
     * \return The position of the i-th one in the outdegrees vector.
     */
    int64_t outdegs_select1(int64_t i) const {return outdegs_ss1.select(i);}

    /**
     * +param i
     * \return The position of the i-th zero in the degrees vector.
     */
    int64_t outdegs_select0(int64_t i) const {return outdegs_ss0.select(i);}

    /**
     * 
     * \return A copy of the outlabels string of the index.
     */
    string get_outlabels() const{
        string S;
        for(int64_t i = 0; i < (LL)outlabels.size(); i++) S += outlabels[i];
        return S;   
    }

    /**
     * \return The size of the outlabels string of the index. 
     */
    int64_t outlabels_size() const {return outlabels.size();}

    /**
    * \param i Index in the outlabels string.
    * \return The character at index i of the outlabels string.
    */ 
    char outlabels_at(int64_t i) const {return outlabels[i];}

    /**
     * \param i Index in the outlabels string.
     * \param c A character.
     * \return Number of times c occurs in outlabels in the half-open range [0,i).
     */
    int64_t outlabels_rank(int64_t i, char c) const {return outlabels.rank(i, c);}

    /**
     * \param i
     * \param c A character.
     * \return The position of the i-th occurrence of character c.
     */
    int64_t outlabels_select(int64_t i, char c) const {return outlabels.select(i, c);}

    /**
     * \return A copy of the C-array of the index.
     */
    vector<int64_t > get_C_array() const{return C;}

    /**
     * \param c
     * \return The C-array at index i, that is, the number of characters in the index with value strictly smaller than c.
     */
    int64_t C_array_at(char c) const {return C[c];}

    /**
     * \return The size of the alphabet of the index.
     */
    int64_t alphabet_size() const {return alphabet.size();}

    /**
     * \return A copy of the alphabet of the index.
     */
    vector<char> get_alphabet() const {return alphabet;}

    /**
     * \param i
     * \return The i-th character of the alphabet of the index.
     */
    char alphabet_at(int64_t i) const {return alphabet[i];}

    /**
     * Serializes the index to the given output stream. The output stream
     * must be in binary mode!
     * 
     * \return The number of bytes written,
     * \param os The output stream
     */
    LL serialize(ostream& os) const{
        LL written = 0;
        written += outlabels.serialize(os);
        written += indegs.serialize(os);
        written += outdegs.serialize(os);

        written += indegs_rs.serialize(os);
        written += indegs_ss1.serialize(os);
        written += indegs_ss0.serialize(os);

        written += outdegs_rs.serialize(os);
        written += outdegs_ss1.serialize(os);
        written += outdegs_ss0.serialize(os);

        // Write C-array
        LL C_array_n_bytes = sizeof(int64_t) * C.size();
        os.write((char*)&C_array_n_bytes, sizeof(C_array_n_bytes));
        os.write((char*)C.data(), C_array_n_bytes);
        written += sizeof(C_array_n_bytes) + C_array_n_bytes;

        // Write alphabet
        LL alphabet_n_bytes = sizeof(char) * alphabet.size();
        os.write((char*)&alphabet_n_bytes, sizeof(alphabet_n_bytes));
        os.write((char*)alphabet.data(), alphabet_n_bytes);
        written += sizeof(alphabet_n_bytes) + alphabet_n_bytes;

        // Write number of nodes
        os.write((char*)&n_nodes, sizeof(n_nodes));
        written += sizeof(n_nodes);

        return written;
    }

    /**
     * Loads the index from the given input stream. The stream must be in binary mode!
     * 
     * \param is The input stream
     */
    void load(istream& is){
        outlabels.load(is);
        indegs.load(is);
        outdegs.load(is);

        indegs_rs.load(is);
        indegs_ss1.load(is);
        indegs_ss0.load(is);

        outdegs_rs.load(is);
        outdegs_ss1.load(is);
        outdegs_ss0.load(is);

        LL C_array_n_bytes;
        is.read((char*)&C_array_n_bytes, sizeof(LL));
        C.resize(C_array_n_bytes / sizeof(int64_t));
        is.read((char*)C.data(), C_array_n_bytes);

        LL alphabet_n_bytes;
        is.read((char*)&alphabet_n_bytes, sizeof(LL));
        alphabet.resize(alphabet_n_bytes / sizeof(char));
        is.read((char*)alphabet.data(), alphabet_n_bytes);

        is.read((char*)&n_nodes, sizeof(n_nodes));

        set_supports();
    }



    WheelerIndex(const WheelerIndex& other){
        assert(&other != this); // What on earth are you trying to do?
        operator=(other);
    }

    WheelerIndex& operator=(const WheelerIndex& other){
        if(&other != this){
            this->outlabels = other.outlabels;

            this->indegs = other.indegs;
            this->outdegs = other.outdegs;
            this->C = other.C;

            this->indegs_rs = other.indegs_rs;
            this->indegs_ss1 = other.indegs_ss1;
            this->indegs_ss0 = other.indegs_ss0;

            this->outdegs_rs = other.outdegs_rs;
            this->outdegs_ss1 = other.outdegs_ss1;
            this->outdegs_ss0 = other.outdegs_ss0;
            
            this->alphabet = other.alphabet;
            this->n_nodes = other.n_nodes;

            set_supports();
            return *this;
        } else return *this; // Assignment to self -> do nothing.
    }

    WheelerIndex() : WheelerIndex(EMPTY_WI_OUTLABELS, EMPTY_WI_INDEGS, EMPTY_WI_OUTDEGS, EMPTY_WI_C) {}

    /**
     * Constructs the index from the given components.
     */
    WheelerIndex(const string& outlabels_string, const sdsl::bit_vector& indegs, const sdsl::bit_vector& outdegs, const vector<int64_t >& C) 
    : indegs(indegs), outdegs(outdegs), C(C){
        sdsl::construct_im(outlabels, outlabels_string.c_str(), 1); // 1: file format is a sequence, not a serialized sdsl object
        sdsl::util::init_support(indegs_rs, &(this->indegs));
        sdsl::util::init_support(indegs_ss1, &(this->indegs));
        sdsl::util::init_support(indegs_ss0, &(this->indegs));
        sdsl::util::init_support(outdegs_rs, &(this->outdegs));
        sdsl::util::init_support(outdegs_ss1, &(this->outdegs));
        sdsl::util::init_support(outdegs_ss0, &(this->outdegs));
        n_nodes = indegs_rs(indegs.size());
        set<char> alphabet_set(outlabels_string.begin(), outlabels_string.end());
        for(char c : alphabet_set) alphabet.push_back(c);
    }

    ~WheelerIndex(){}

    /**
     *  \param edge The Wheeler rank of an edge.
     *  \return The Wheeler rank of the node at the destination of the edge.
     */
    int64_t edge_destination(int64_t edge) const{
        assert(edge >= 0 && edge < (LL)outlabels.size());
        return indegs_rank1(indegs_select0(edge+1)) - 1;
    }

    /**
     *  \param edge The Wheeler rank of an edge.
     *  \return The Wheeler rank of the node at the source of the edge.
     */
    int64_t edge_source(int64_t edge) const{
        assert(edge >= 0 && edge < outlabels.size());
        LL c = edge_label(edge);
        LL outlabels_idx = outlabels_select(edge - C[c] + 1, c);
        return outdegs_rank1(outdegs_select0(outlabels_idx+1)) - 1;
    }

    /**
     * \param edge The Wheeler rank of an edge.
     * \return The label of the edge.
     * \note Has time complexity O(alphabet size).
     */
    char edge_label(int64_t edge) const{
        // Todo: maybe use a lookup table for constant time commplexity
        char c = 0;
        for(char a : alphabet) if(C[a] <= edge) c = a;
        return c;
    }

    /**
     * \param l
     * \param r
     * \param c
     * \return Given a node range [l,r] and a character c, returns the range of outedges from [l,r] with label c. If no such edge exists, returns (1,0).
     */
    pair<int64_t, int64_t> nodes_to_outedges(int64_t l, int64_t r, char c) const{
        int64_t start = outdegs_select1(l+1) - l; // Start of the outlabels in the search range
        int64_t end;
        if(r == n_nodes-1) end = outlabels_size()-1;
        else end = outdegs_select1(r+1+1) - (r + 1) - 1; // Inclusive end of the outlabels in the search range
        if(end < start) return {1,0}; // No outlabels in the range

        // [start,end] is the range of outgoing labels in outlabels
        int64_t edge_l = outlabels_rank(start, c); // number of c-edges before the range
        int64_t edge_r = outlabels_rank(end+1, c); // number of c-edges up to the end of the range
        int64_t num_c_edges = edge_r - edge_l;

        if(num_c_edges == 0) return {1,0}; // No outgoing c

        return {C[c] + edge_l, C[c] + edge_r - 1};
    }

    /**
     * \param l
     * \param r
     * \return Given an edge Wheeler range [l,r], returns the range of nodes at the destinations of the edges.
     */
    pair<int64_t, int64_t> inedges_to_nodes(int64_t l, int64_t r) const{
        if(l > r) return {1,0};
        return {edge_destination(l), edge_destination(r)};
    }

    /**
     * \param l
     * \param r
     * \param c
     * \return The range of nodes that are reached by following an edge labeled with c from the Wheeler range of nodes [l,r].
     *         If the resulting range is empty, returns (1,0).
     */
    pair<int64_t, int64_t> node_range_follow_edges(int64_t l, int64_t r, int64_t c) const{
        std::tie(l,r) = nodes_to_outedges(l,r,c);
        if(l > r) return {1,0};
        return inedges_to_nodes(l,r);
    }

    /**
     * \param S Search string.
     * \param len Length of the search string.
     * \return The range of nodes with incoming path label S, or (1,0) if no such node exists.
     */
    pair<int64_t, int64_t> search(const char* S, int64_t len) const{
        // [l,r] = current colexicographic range (zero-indexed)
        int64_t l, r;
        std::tie(l,r) = full_node_range();

        for(int64_t i = 0; i < len; i++){
            std::tie(l,r) = node_range_follow_edges(l,r,S[i]);
            if(r < l) return {1,0};
        }
        return {l,r};
    }

    /**
     * \param S Search string.
     * \return The range of nodes with incoming path label S, or (1,0) if no such node exists.
     */
    pair<int64_t, int64_t> search(const string& S) const {
        return search(S.c_str(), S.size());
    }

    /**
     * \param node
     * \return A range [l,r] of indices in outlabels of the outgoing edges from node. If the node has no outgoing labels, returns (1,0).
     */
    pair<int64_t, int64_t> outedge_range(int64_t node) const{
        if(node == -1) return {1,0};

        // Check if the node has an outedge with c. If no, return fail. If yes, we need to find the Wheeler rank 
        // of that outedge. Then we use the indegree vector to find the destination of that edge
        // See also:
        // Alanko, J., Gagie, T., Navarro, G., & Benkner, L. S. (2018). Tunneling on Wheeler Graphs.
        // arXiv preprint arXiv:1811.02457.
        
        int64_t outdegs_idx = outdegs_ss1(node+1);
        int64_t outdegree = 0;
        while(outdegs_idx + outdegree + 1 <= (LL)outdegs.size()-1 && outdegs[outdegs_idx + outdegree + 1] == 0)
            outdegree++;
        if(outdegree == 0) return {1,0};

        int64_t start = outdegs_idx - node;
        int64_t end = start + outdegree - 1;
        return {start,end};
    }

    /**
     * \param node
     * \return A range [l,r] of the Wheeler ranks of the incoming edges to a node. If the node has no incoming labels, returns (1,0).
     */
    pair<int64_t, int64_t> inedge_range(int64_t node) const{
        if(node == -1) return {1,0};
        LL indegs_idx = indegs_select1(node + 1);
        LL indegree = this->indegree(node); // Todo optimize: can save one select by reusing indegs_idx
        if(indegree == 0) return {1,0};
        LL start = indegs_rank0(indegs_idx);
        return {start, start + indegree - 1};
    }
    
    /**
     * \param node
     * \param answer A buffer for the output.
     * \return Stores the outgoing edge labels from node into the given character array. The caller should initialize the array to have at least one char of space for each character in the alphabet. Does not add a *null* at the end of answer.
     */
    int64_t node_outlabels(int64_t node, char* answer) const{
        if(node == -1) return 0;

        int64_t start, end;
        std::tie(start,end) = outedge_range(node);
        if(end < start) return 0;

        int64_t count = 0;
        for(int64_t i = 0; i < end-start+1; i++){
            answer[i] = outlabels[start + i];
            count++;
        }
        return count;
    }

    /**
     * \param node
     * \return The edge labels of the outgoing edges from the node.
     */
    string node_outlabels(int64_t node) const{
        char answer[alphabet.size()];
        int64_t len = node_outlabels(node, answer);
        string S(answer, answer + len);
        return S;
    }

    /**
     * \param node
     * \return The number of incoming edges to node.
     * \note Takes time proportional to the indegree.
     */
    int64_t indegree(int64_t node) const{
        // todo: unit test
        if(node == -1) return -1;

        int64_t idx = indegs_ss1(node+1) + 1;
        int64_t ans = 0;
        while(idx < (LL)indegs.size() && indegs[idx] == 0){
            idx++;
            ans++;
        }
        return ans;
    }

    /**
     * \param node
     * \return The number of outgoing edges from node.
     * \note Takes time proportional to the outdegree.
     */
    int64_t outdegree(int64_t node) const{
        // todo: unit test
        if(node == -1) return -1;

        int64_t idx = outdegs_ss1(node+1) + 1;
        int64_t ans = 0;
        while(idx < (LL)outdegs.size() && outdegs[idx] == 0){
            idx++;
            ans++;
        }
        return ans;
    }

    /**
     * \param node
     * \return The edge label of the incoming characters to the node. Note that since the graph is Wheeler, all incoming edges to a node will have the same label.
     */
    char incoming_character(int64_t node) const{ // todo: test
        assert(node != -1);

        int64_t pos_in_indegs = indegs_ss1(node+1);
        if(pos_in_indegs == (LL)indegs.size()-1 || indegs[pos_in_indegs+1] == 1){
            // Indegree zero.
            return '\0';
        } 
        int64_t edge_rank = pos_in_indegs - node; // zero-based rank
        return edge_label(edge_rank);
    }

    /**
     * \param idx Edge index in outlabels.
     * \return The Wheeler rank of the edge.
     */
    LL outedge_index_to_wheeler_rank(int64_t idx) const{
        char c = outlabels_at(idx);
        return C_array_at(c) + outlabels_rank(idx,c);
    }

    string to_string() const{ // For debugging
        string bwt = get_outlabels();
        string bwt_spaced;

        LL next = 0;
        for(LL i = 0; i < outdegs.size(); i++){
            if(outdegs[i] == 1) bwt_spaced += " ";
            else{
                bwt_spaced += bwt[next];
                next++;
            }
        }
        
        stringstream ss;
        ss << bwt_spaced << "\n" << outdegs << "\n" << indegs << "\n";
        ss << "C = ";
        for(LL i = 0; i < C.size(); i++){
            if(i == 0 || C[i-1] != C[i]) ss<< C[i] << " "; 
        }

        return ss.str();
    }

    template<typename T>
    friend bool operator==(const WheelerIndex<T>& wi1, const WheelerIndex<T>& wi2);
};

template<typename bitvector_t>
bool operator==(const WheelerIndex<bitvector_t>& wi1, const WheelerIndex<bitvector_t>& wi2){
    // Returns true if the data is the same and the rank and
    // select structures give identical answers to aint64_t queries

    bool same = true;

    if(wi1.indegs_size() != wi2.indegs_size()) return false;
    for(int64_t i = 0; i < wi1.indegs_size(); i++)
        same &= wi1.indegs_at(i) == wi2.indegs_at(i);

    if(wi1.outdegs_size() != wi2.outdegs_size()) return false;
    for(int64_t i = 0; i < wi1.outdegs_size(); i++)
        same &= wi1.outdegs_at(i) == wi2.outdegs_at(i);

    if(wi1.get_outlabels() != wi2.get_outlabels()) return false;
    if(wi1.C != wi2.C) return false;
    if(wi1.number_of_nodes() != wi2.number_of_nodes()) return false;

    for(int64_t i = 0; i <= (LL)wi1.indegs.size(); i++){
        if(wi1.indegs_rs(i) != wi2.indegs_rs(i)) return false;
    }

    for(int64_t i = 1; i <= wi1.number_of_nodes(); i++){
        if(wi1.indegs_ss1(i) != wi2.indegs_ss1(i)) return false;
    }

    for(int64_t i = 1; i <= (LL)wi1.indegs.size() - wi1.number_of_nodes(); i++){
        if(wi1.indegs_ss0(i) != wi2.indegs_ss0(i)) return false;
    }

    for(int64_t i = 1; i <= wi1.number_of_nodes(); i++){
        if(wi1.outdegs_rs(i) != wi2.outdegs_rs(i)) return false;
    }

    for(int64_t i = 1; i <= wi1.number_of_nodes(); i++){
        if(wi1.outdegs_ss1(i) != wi2.outdegs_ss1(i)) return false;
    }

    for(int64_t i = 1; i <= (LL)wi1.outdegs.size() - wi1.number_of_nodes(); i++){
        if(wi1.outdegs_ss0(i) != wi2.outdegs_ss0(i)) return false;
    }

    if(wi1.alphabet != wi2.alphabet) return false;

    return true;
}

} // End of namespace