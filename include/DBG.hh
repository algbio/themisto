#pragma once

#include "libwheeler/BOSS.hh"

class DBG{

public:

BOSS<sdsl::bit_vector>* boss;
vector<bool> is_dummy;

    struct Node{
        LL id;
    };

    struct Edge{
        LL source;
        LL dest;
        char label;
    };

    class all_nodes_generator;
    class outedge_iterator;
    class inedge_iterator;

    DBG(BOSS<sdsl::bit_vector>* boss) : boss(boss), is_dummy(boss->get_dummy_node_marks()){}

    all_nodes_generator all_nodes();
    outedge_iterator outedges(Node v);
    inedge_iterator inedges(Node v);

};


class DBG::all_nodes_generator{

public:

    struct iterator{

        int idx;
        BOSS<sdsl::bit_vector>* boss;
        vector<bool>* is_dummy;

        iterator(int node_idx, BOSS<sdsl::bit_vector>* boss, vector<bool>* is_dummy) : idx(node_idx), boss(boss), is_dummy(is_dummy) {
            // Rewind to first non-dummy
            while(idx < boss->number_of_nodes() && is_dummy->at(idx)) idx++;
        }

        iterator operator++(){
            idx++;
            while(idx < boss->number_of_nodes() && is_dummy->at(idx)) idx++;
            return *this;
        }

        Node operator*(){
            return {idx};
        }

        bool operator!=(iterator& other){
            return this->idx != other.idx;
        }
    };

    BOSS<sdsl::bit_vector>* boss;
    vector<bool>* is_dummy;
    all_nodes_generator(BOSS<sdsl::bit_vector>* boss, vector<bool>* is_dummy) : boss(boss), is_dummy(is_dummy){}

    iterator begin(){return iterator(0, boss, is_dummy);}
    iterator end(){return iterator(boss->number_of_nodes(), boss, is_dummy);}

};

DBG::all_nodes_generator DBG::all_nodes(){
    return all_nodes_generator(boss, &is_dummy);
}


/*

for(DBG::Node node : DBG){
    for(DBG::Edge edge : DBG.outedges(node)){
        DBG::Node destination = edge.destination();
    }
}

*/