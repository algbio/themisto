#pragma once

#include <unordered_map>
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"
#include "backward_traversal.hh"

using namespace std;
using namespace sbwt;

// A class that offers an interface to the de Bruijn graph.
// Implemented internally using the SBWT class. The SBWT is a bit tricky
// because it has the dummy nodes. This class deals with the dummies so 
// that the caller does not have to care about them.

/*

Usage:

for(DBG::Node node : DBG){
    for(DBG::Edge edge : DBG.outedges(node)){
        DBG::Node destination = edge.dest
    }
}

*/

class DBG{

private:

    // No copying because of pointer business
    DBG(DBG const& other) = delete;
    DBG(DBG&) = default;

    const plain_matrix_sbwt_t* SBWT; // Non-owning pointer
    SBWT_backward_traversal_support* backward_support; // Owning pointer

    // Rank support to the internal dummy bit vector of backward_support. This is
    // a little unsafe I know, but it's probably ok because this object controls the memory
    // lifetime of the backward support.
    sdsl::rank_support_v5<> dummy_node_rs; 

public:

    struct Node{
        int64_t id; // From 0 to number of subsets in the SBWT.

        Node(int64_t id) : id(id) {}

        bool operator==(const Node& other) const{
            return this->id == other.id;
        }
    };

    struct Edge{ 
        int64_t source;
        int64_t dest;
        char label;
    };

    class all_nodes_generator;
    class outedge_generator;
    class inedge_generator;

    DBG(){}
    DBG(const plain_matrix_sbwt_t* SBWT) : SBWT(SBWT){
        backward_support = new SBWT_backward_traversal_support(SBWT);
        sdsl::util::init_support(dummy_node_rs, &backward_support->get_dummy_marks());
    }

    all_nodes_generator all_nodes() const; // Return a generator for a range-for loop
    outedge_generator outedges(Node v) const; // Return a generator for a range-for loop
    inedge_generator inedges(Node v) const; // Return a generator for a range-for loop

    Node locate(const string& kmer) const{ // Returns a node with id -1 if the k-mer does not does not exist
        return {SBWT->search(kmer)};
    }

    string get_node_label(const Node& v) const{
        assert(v.id != -1);
        return backward_support->get_node_label(v.id);
    }

    int64_t indegree(const Node& v) const{
        assert(v.id != -1);
        int64_t in_neighbors[4]; 
        int64_t indeg;
        backward_support->list_DBG_in_neighbors(v.id, in_neighbors, indeg);
        return indeg;
    }

    int64_t outdegree(const Node& v) const{
        assert(v.id != -1);
        int64_t node = v.id;
        while(SBWT->get_streaming_support()[node] == 0) node--; // Walk back to the start of the suffix group
        return SBWT->get_subset_rank_structure().A_bits[node] 
             + SBWT->get_subset_rank_structure().C_bits[node]
             + SBWT->get_subset_rank_structure().G_bits[node]
             + SBWT->get_subset_rank_structure().T_bits[node];
    }

    // If the indegree of v is not 1, throws an error.
    // Otherwise, returns the node at the start of the incoming edge to v
    Node pred(const Node& v) const{
        if(indegree(v) != 1) 
            throw std::invalid_argument("Tried to get the predecessor of a node with indegree " + to_string(indegree(v)));
        else{
            return {backward_support->backward_step(v.id)};
        }
    }

    // If the outdegree of v is not 1, throws an error.
    // Otherwise, returns the node at the end of the outgoing edge to v
    Node succ(const Node& v) const{
        int64_t node = v.id;
        while(SBWT->get_streaming_support()[node] == 0) node--; // Walk back to the start of the suffix group

        if(SBWT->get_subset_rank_structure().A_bits[node]) return {SBWT->forward(node, 'A')};
        else if(SBWT->get_subset_rank_structure().C_bits[node]) return {SBWT->forward(node, 'C')};
        else if(SBWT->get_subset_rank_structure().G_bits[node]) return {SBWT->forward(node, 'G')};
        else if(SBWT->get_subset_rank_structure().T_bits[node]) return {SBWT->forward(node, 'T')};
        else throw std::invalid_argument("Tried to get the successor of a node with outdegree " + to_string(outdegree(v)));
    }

    char incoming_character(const Node& v) const{
        return backward_support->get_incoming_character(v.id);
    }


    int64_t number_of_kmers() const{
        return SBWT->number_of_kmers();
    }

    int64_t get_k() const{
        return SBWT->get_k();
    }

    // If v is the i-th k-mer in colexicographic order (0-based), returns i.
    int64_t kmer_rank(Node v) const{
        if(v.id == -1) return -1;
        return v.id - dummy_node_rs.rank(v.id); // Subtract dummy nodes
    }

    int64_t number_of_sets_in_sbwt() const {
        return SBWT->number_of_subsets();
    }

    // Low level function for people who know what they are doing
    int64_t is_dummy_colex_position(int64_t colex) const {
        return backward_support->get_dummy_marks()[colex];
    }

    ~DBG(){
        delete backward_support;
    }

};


class DBG::all_nodes_generator{

public:

    struct end_iterator{}; // Dummy end iterator
    struct iterator{

        int64_t idx;
        const plain_matrix_sbwt_t* SBWT;
        const sdsl::bit_vector* dummy_marks;

        iterator(int64_t node_idx, const plain_matrix_sbwt_t* SBWT, const sdsl::bit_vector* dummy_marks) : idx(node_idx), SBWT(SBWT), dummy_marks(dummy_marks){
            // Rewind to first non-dummy
            while(idx < SBWT->number_of_subsets() && (*dummy_marks)[idx]) idx++;
        }

        iterator operator++(){
            idx++;
            while(idx < SBWT->number_of_subsets() && (*dummy_marks)[idx]) idx++;
            return *this;
        }

        Node operator*(){
            return Node(idx);
        }

        bool operator!=(const end_iterator& other){
            (void)other; // This is just a dummy
            return this->idx < SBWT->number_of_subsets();
        }
    };

    const plain_matrix_sbwt_t* SBWT;
    const sdsl::bit_vector* dummy_marks;
    all_nodes_generator(const plain_matrix_sbwt_t* SBWT, const sdsl::bit_vector* dummy_marks) : SBWT(SBWT), dummy_marks(dummy_marks){}

    iterator begin(){return iterator(0, SBWT, dummy_marks);}
    end_iterator end(){return end_iterator();}

};


class DBG::outedge_generator{

public:

    struct end_iterator{}; // Dummy end iterator
    struct iterator{

        int64_t source_node; // For reporting the source node in the edge
        int64_t suffix_group_start;
        int64_t sbwt_row;
        const plain_matrix_sbwt_t* SBWT;

        iterator(int64_t node_idx, const plain_matrix_sbwt_t* SBWT) : source_node(node_idx), suffix_group_start(node_idx), sbwt_row(-1), SBWT(SBWT){
            while(SBWT->get_streaming_support()[suffix_group_start] == 0) suffix_group_start--; // Walk back to the start of the suffix group
            operator++(); // Rewind to the first outedge
        }

        bool matrixboss_access(int64_t row, int64_t col) const{
            assert(row >= 0 && row < 4);
            if(row == 0) return SBWT->get_subset_rank_structure().A_bits[col];
            if(row == 1) return SBWT->get_subset_rank_structure().C_bits[col];
            if(row == 2) return SBWT->get_subset_rank_structure().G_bits[col];
            if(row == 3) return SBWT->get_subset_rank_structure().T_bits[col];
            return false; // SHould not happen
        }

        iterator operator++(){
            sbwt_row++;
            while(sbwt_row < 4 && matrixboss_access(sbwt_row, suffix_group_start) == 0) sbwt_row++;
            return *this;
        }

        Edge operator*(){
            char label = sbwt::char_idx_to_DNA(sbwt_row);
            int64_t dest = SBWT->forward(suffix_group_start, label);
            return {.source = source_node, .dest = dest, .label = label};
        }

        bool operator!=(const end_iterator& other){
            (void)other; // This is just a dummy
            return sbwt_row < 4;
        }
    };

    Node v;
    const plain_matrix_sbwt_t* SBWT;

    outedge_generator(Node v, const plain_matrix_sbwt_t* SBWT) : v(v), SBWT(SBWT){}

    iterator begin(){return iterator(v.id, SBWT);}
    end_iterator end(){return end_iterator();}

};


class DBG::inedge_generator{

public:

    struct end_iterator{}; // Dummy end iterator
    struct iterator{

        int64_t node_idx;
        const plain_matrix_sbwt_t* SBWT;
        const SBWT_backward_traversal_support* backward_support;

        int64_t in_neighbors[4];
        int64_t indegree = 0;
        int64_t in_neighbors_idx = 0;
        char incoming_char = 0;

        iterator(int64_t node_idx, const plain_matrix_sbwt_t* SBWT, const SBWT_backward_traversal_support* backward_support) : node_idx(node_idx), SBWT(SBWT), backward_support(backward_support) {
            backward_support->list_DBG_in_neighbors(node_idx, in_neighbors, indegree);
            incoming_char = backward_support->get_incoming_character(node_idx);
        }

        iterator operator++(){
            in_neighbors_idx++;
            return *this;
        }

        Edge operator*(){
            int64_t dest = node_idx;
            int64_t source = in_neighbors[in_neighbors_idx];
            return {.source = source, .dest = dest, .label = incoming_char};
        }

        bool operator!=(const end_iterator& other){
            (void)other; // This is just a dummy
            return in_neighbors_idx < indegree;
        }
    };

    Node v;
    const plain_matrix_sbwt_t* SBWT;
    const SBWT_backward_traversal_support* backward_support;

    inedge_generator(Node v, const plain_matrix_sbwt_t* SBWT, const SBWT_backward_traversal_support* backward_support) : v(v), SBWT(SBWT), backward_support(backward_support) {}

    iterator begin(){return iterator(v.id, SBWT, backward_support);}
    end_iterator end(){return end_iterator();}

};

// For standard library hash functions.
template<>
struct std::hash<DBG::Node>{
    std::size_t operator()(const DBG::Node& X) const {
        return X.id;
    }
};

