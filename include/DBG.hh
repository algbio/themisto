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
    class outedge_generator;
    class inedge_generator;

    DBG(BOSS<sdsl::bit_vector>* boss) : boss(boss), is_dummy(boss->get_dummy_node_marks()){}

    all_nodes_generator all_nodes();
    outedge_generator outedges(Node v);
    inedge_generator inedges(Node v);

};


class DBG::all_nodes_generator{

public:

    struct iterator{

        LL idx;
        BOSS<sdsl::bit_vector>* boss;
        vector<bool>* is_dummy;

        iterator(LL node_idx, BOSS<sdsl::bit_vector>* boss, vector<bool>* is_dummy) : idx(node_idx), boss(boss), is_dummy(is_dummy) {
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


class DBG::outedge_generator{

public:

    struct iterator{

        LL node_idx;
        LL outlabels_start; // In outlabels of boss
        LL outlabels_offset; // In outlabels of boss
        BOSS<sdsl::bit_vector>* boss;
        vector<bool>* is_dummy;

        iterator(LL node_idx, LL edge_offset, BOSS<sdsl::bit_vector>* boss, vector<bool>* is_dummy) : node_idx(node_idx), outlabels_offset(edge_offset), boss(boss), is_dummy(is_dummy) {
            outlabels_start = boss->outdegs_rank0(boss->outdegs_select1(node_idx+1));
        }

        iterator operator++(){
            outlabels_offset++;
            return *this;
        }

        Edge operator*(){
            LL source = node_idx;
            char label = boss->outlabels_at(outlabels_start + outlabels_offset);
            LL dest = boss->edge_destination(boss->outedge_index_to_wheeler_rank(outlabels_start + outlabels_offset));
            return {.source = source, .dest = dest, .label = label};
        }

        bool operator!=(iterator& other){
            return this->outlabels_offset != other.outlabels_offset;
        }
    };

    Node v;
    BOSS<sdsl::bit_vector>* boss;
    vector<bool>* is_dummy;

    outedge_generator(Node v, BOSS<sdsl::bit_vector>* boss, vector<bool>* is_dummy) : v(v), boss(boss), is_dummy(is_dummy){}

    iterator begin(){return iterator(v.id, 0, boss, is_dummy);}
    iterator end(){return iterator(v.id, boss->outdegree(v.id), boss, is_dummy);}

};


DBG::all_nodes_generator DBG::all_nodes(){
    return all_nodes_generator(boss, &is_dummy);
}

DBG::outedge_generator DBG::outedges(Node v){
    return outedge_generator(v, boss, &is_dummy);
}


/*

for(DBG::Node node : DBG){
    for(DBG::Edge edge : DBG.outedges(node)){
        DBG::Node destination = edge.destination();
    }
}

*/