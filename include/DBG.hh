#pragma once

#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"
#include "backward_traversal.hh"

using namespace std;
using namespace sbwt;

// A class that offers an interface to the de Bruijn graph.
// Implemented internally using the SBWT class. The SBWT is a bit tricky
// because it has the dummy nodes. This class deals with the dummies so 
// that the caller does not have to care about them.
class DBG{

private:

    const plain_matrix_sbwt_t* SBWT; // Non-owning pointer
    SBWT_backward_traversal_support backward_support;

public:

    struct Node{
        int64_t id;

        Node(int64_t id) : id(id) {}
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
    DBG(const plain_matrix_sbwt_t* SBWT) : SBWT(SBWT), backward_support(SBWT){}

    all_nodes_generator all_nodes(); // Return a generator for a range-for loop
    outedge_generator outedges(Node v); // Return a generator for a range-for loop
    inedge_generator inedges(Node v); // Return a generator for a range-for loop

    Node locate(const string& kmer){ // Returns a node with id -1 if the k-mer does not does not exist
        return {SBWT->search(kmer)};
    }

    string get_node_label(const Node& v){
        assert(v.id != -1);
        return backward_support.get_node_label(v.id);
    }

    int64_t indegree(const Node& v){
        assert(v.id != -1);
        int64_t in_neighbors[4]; 
        int64_t indeg;
        backward_support.list_DBG_in_neighbors(v.id, in_neighbors, indeg);
        return indeg;
    }

    int64_t outdegree(const Node& v){
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
    Node pred(const Node& v){
        if(indegree(v) != 1) 
            throw std::invalid_argument("Tried to get the predecessor of a node with indegree " + to_string(indegree(v)));
        else{
            return {backward_support.backward_step(v.id)};
        }
    }

    // If the outdegree of v is not 1, throws an error.
    // Otherwise, returns the node at the end of the outgoing edge to v
    Node succ(const Node& v){
        int64_t node = v.id;
        while(SBWT->get_streaming_support()[node] == 0) node--; // Walk back to the start of the suffix group

        if(SBWT->get_subset_rank_structure().A_bits[node]) return {SBWT->forward(node, 'A')};
        else if(SBWT->get_subset_rank_structure().C_bits[node]) return {SBWT->forward(node, 'C')};
        else if(SBWT->get_subset_rank_structure().G_bits[node]) return {SBWT->forward(node, 'G')};
        else if(SBWT->get_subset_rank_structure().T_bits[node]) return {SBWT->forward(node, 'T')};
        else throw std::invalid_argument("Tried to get the predecessor of a node with indegree " + to_string(indegree(v)));
    }

    char incoming_character(const Node& v){
        return backward_support.get_incoming_character(v.id);
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
            return {.id = idx};
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

        int64_t node_idx;
        int64_t sbwt_row;
        const plain_matrix_sbwt_t* SBWT;

        iterator(int64_t node_idx, const plain_matrix_sbwt_t* SBWT) : node_idx(node_idx), sbwt_row(-1), SBWT(SBWT){
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
            while(sbwt_row < 4 && matrixboss_access(sbwt_row, node_idx) == 0) sbwt_row++;
            return *this;
        }

        Edge operator*(){
            int64_t source = node_idx;
            char label = char_idx_to_DNA(sbwt_row);
            int64_t dest = SBWT->forward(source, label);
            return {.source = source, .dest = dest, .label = label};
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



DBG::all_nodes_generator DBG::all_nodes(){
    return all_nodes_generator(SBWT, &(backward_support.get_dummy_marks()));
}

DBG::outedge_generator DBG::outedges(Node v){
    return outedge_generator(v, SBWT);
}

DBG::inedge_generator DBG::inedges(Node v){
    return inedge_generator(v, SBWT, &backward_support);
}

// For standard library hash functions.
struct std::hash<DBG::Node>{
    std::size_t operator()(const DBG::Node& X) const {
        return X.id;
    }
};


/*

Usage:

for(DBG::Node node : DBG){
    for(DBG::Edge edge : DBG.outedges(node)){
        DBG::Node destination = edge.dest
    }
}

*/