#pragma once

#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"
#include "backward_traversal.hh"

using namespace std;
using namespace sbwt;

// A class that offers an interface to the de Bruijn graph.
// Implemented internally using the BOSS class. The BOSS is a bit tricky
// because it has the dummy nodes. This class deals with the dummies so 
// that the caller does not have to care about them.
class DBG{

private:

    plain_matrix_sbwt_t* SBWT; // Non-owning pointer
    SBWT_backward_traversal_support backward_support;
    vector<bool> is_dummy;

public:

    struct Node{
        int64_t id;
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
    DBG(plain_matrix_sbwt_t* SBWT) : SBWT(SBWT), backward_support(SBWT), is_dummy(SBWT->compute_dummy_node_marks()){}

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
        if(v.id == 0) return 0; // Root
        int64_t u = backward_support...
        pair<LL,LL> R = boss->inedge_range(v.id);
        if(R.first == R.second){ // One predecessor. Check whether it's a dummy
            int64_t source = boss->edge_source(R.first);
            if(is_dummy.at(source)) return 0;
        } 
        return R.second - R.first + 1;
    }

    int64_t outdegree(const Node& v){
        assert(v.id != -1);
        return boss->outdegree(v.id);
    }

    // If the indegree of v is not 1, throws an error.
    // Otherwise, returns the node at the start of the incoming edge to v
    Node pred(const Node& v){
        if(indegree(v) != 1) 
            throw std::invalid_argument("Tried to get the predecessor of a node with indegree " + to_string(indegree(v)));
        else{
            return {boss->edge_source(boss->inedge_range(v.id).first)}; // predecessor
        }
    }

    // If the outdegree of v is not 1, throws an error.
    // Otherwise, returns the node at the end of the outgoing edge to v
    Node succ(const Node& v){
        if(outdegree(v) != 1) 
            throw std::invalid_argument("Tried to get the predecessor of a node with indegree " + to_string(indegree(v)));
        else{
            return {boss->walk(v.id, boss->outlabels_at(boss->outlabel_range(v.id).first))}; // successor
        }
    }    

};


class DBG::all_nodes_generator{

public:

    struct end_iterator{}; // Dummy end iterator
    struct iterator{

        int64_t idx;
        BOSS<sdsl::bit_vector>* boss;
        vector<bool>* is_dummy;

        iterator(int64_t node_idx, BOSS<sdsl::bit_vector>* boss, vector<bool>* is_dummy) : idx(node_idx), boss(boss), is_dummy(is_dummy) {
            // Rewind to first non-dummy
            while(idx < boss->number_of_nodes() && is_dummy->at(idx)) idx++;
        }

        iterator operator++(){
            idx++;
            while(idx < boss->number_of_nodes() && is_dummy->at(idx)) idx++;
            return *this;
        }

        Node operator*(){
            return {.id = idx};
        }

        bool operator!=(const end_iterator& other){
            (void)other; // This is just a dummy
            return this->idx < boss->number_of_nodes();
        }
    };

    BOSS<sdsl::bit_vector>* boss;
    vector<bool>* is_dummy;
    all_nodes_generator(BOSS<sdsl::bit_vector>* boss, vector<bool>* is_dummy) : boss(boss), is_dummy(is_dummy){}

    iterator begin(){return iterator(0, boss, is_dummy);}
    end_iterator end(){return end_iterator();}

};


class DBG::outedge_generator{

public:

    struct end_iterator{}; // Dummy end iterator
    struct iterator{

        int64_t node_idx;
        int64_t outlabels_start; // In outlabels of boss
        int64_t outlabels_offset; // In outlabels of boss
        int64_t outdegree;
        BOSS<sdsl::bit_vector>* boss;

        iterator(int64_t node_idx, int64_t edge_offset, BOSS<sdsl::bit_vector>* boss) : node_idx(node_idx), outlabels_offset(edge_offset), boss(boss){
            outlabels_start = boss->outdegs_rank0(boss->outdegs_select1(node_idx+1));
            outdegree = boss->outdegree(node_idx);
        }

        iterator operator++(){
            outlabels_offset++;
            return *this;
        }

        Edge operator*(){
            int64_t source = node_idx;
            char label = boss->outlabels_at(outlabels_start + outlabels_offset);
            int64_t dest = boss->edge_destination(boss->outedge_index_to_wheeler_rank(outlabels_start + outlabels_offset));
            return {.source = source, .dest = dest, .label = label};
        }

        bool operator!=(const end_iterator& other){
            (void)other; // This is just a dummy
            return outlabels_offset < outdegree;
        }
    };

    Node v;
    BOSS<sdsl::bit_vector>* boss;

    outedge_generator(Node v, BOSS<sdsl::bit_vector>* boss) : v(v), boss(boss){}

    iterator begin(){return iterator(v.id, 0, boss);}
    end_iterator end(){return end_iterator();}

};


class DBG::inedge_generator{

// This does a lot of redundant work and could be optimized

public:

    struct end_iterator{}; // Dummy end iterator
    struct iterator{

        int64_t node_idx;
        int64_t inlabels_start; // Wheeler rank in boss
        int64_t inlabels_offset;
        int64_t indegree;
        char incoming_char;
        BOSS<sdsl::bit_vector>* boss;
        vector<bool>* is_dummy;

        iterator(int64_t node_idx, int64_t edge_offset, BOSS<sdsl::bit_vector>* boss, vector<bool>* is_dummy) : node_idx(node_idx), inlabels_offset(edge_offset), boss(boss), is_dummy(is_dummy) {
            inlabels_start = boss->indegs_rank0(boss->indegs_select1(node_idx+1));
            indegree = boss->indegree(node_idx);
            incoming_char = boss->incoming_character(node_idx);
            if(indegree == 1){
                // Check if we are preceeded by a dummy. If yes, there are no DBG in-edges
                int64_t source = boss->edge_source(inlabels_start);
                if(is_dummy->at(source)) inlabels_offset++; // Should match the end iterator now
            }
        }

        iterator operator++(){
            inlabels_offset++;
            return *this;
        }

        Edge operator*(){
            int64_t dest = node_idx;
            int64_t source = boss->edge_source(inlabels_start + inlabels_offset);
            return {.source = source, .dest = dest, .label = incoming_char};
        }

        bool operator!=(const end_iterator& other){
            (void)other; // This is just a dummy
            return inlabels_offset < indegree;
        }
    };

    Node v;
    BOSS<sdsl::bit_vector>* boss;
    vector<bool>* is_dummy;

    inedge_generator(Node v, BOSS<sdsl::bit_vector>* boss, vector<bool>* is_dummy) : v(v), boss(boss), is_dummy(is_dummy){}

    iterator begin(){return iterator(v.id, 0, boss, is_dummy);}
    end_iterator end(){return end_iterator();}

};



DBG::all_nodes_generator DBG::all_nodes(){
    return all_nodes_generator(boss, &is_dummy);
}

DBG::outedge_generator DBG::outedges(Node v){
    return outedge_generator(v, boss);
}

DBG::inedge_generator DBG::inedges(Node v){
    return inedge_generator(v, boss, &is_dummy);
}

/*

Usage:

for(DBG::Node node : DBG){
    for(DBG::Edge edge : DBG.outedges(node)){
        DBG::Node destination = edge.dest
    }
}

*/