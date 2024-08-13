#include <iostream>
#include "DBG.hh"
#include "globals.hh"
#include "coloring/Coloring.hh"

namespace new_extract_unitigs_internals {

bool is_first_kmer_of_unitig(const DBG& dbg, const DBG::Node& node) {
    int64_t indeg = dbg.indegree(node);
    if(indeg != 1){
        return true;
    } else {
        DBG::Node u = dbg.pred(node);
        return dbg.outdegree(u) > 1;
    } }

// Returns the sequence of nodes and the label of the unitig
pair<vector<DBG::Node>, vector<char>> walk_unitig_from(const DBG& dbg, DBG::Node v) {
    DBG::Node v0 = v;
    vector<DBG::Node> nodes;
    nodes.push_back(v);
    
    string first_kmer_label = dbg.get_node_label(v);
    vector<char> label(first_kmer_label.begin(), first_kmer_label.end());

    while(dbg.outdegree(v) == 1) {
        DBG::Edge e = *(dbg.outedges(v).begin());

        char c = e.label;
        v = e.dest; // Follow edge
        if(v != v0 && dbg.indegree(v) == 1) {
            label.push_back(c);
            nodes.push_back(v);
        } else { break; }
    }

    return {nodes, label};
}

static const char* UNITIG_ID_EQ = "unitig_id=";
static const int64_t UNITIG_ID_EQ_LEN = 10;

static const char* COLOR_ID_EQ = "color_id=";
static const int64_t COLOR_ID_EQ_LEN = 9;

void write_unitig(
    const char* unitig_id_chars, int64_t unitig_id_length, 
    const char* unitig_label, int64_t unitig_label_length, 
    const char* color_set_id, int64_t color_set_id_len, 
    ParallelOutputWriter& unitigs_out){

    vector<char> buf; // ASCII data. Needs to be written in one go to avoid interleaved writing with other threads!

    buf.push_back('>');
    buf.insert(buf.end(), UNITIG_ID_EQ, UNITIG_ID_EQ + UNITIG_ID_EQ_LEN);
    buf.insert(buf.end(), unitig_id_chars, unitig_id_chars + unitig_id_length); // Push the unitig id
    buf.push_back(' ');
    buf.insert(buf.end(), COLOR_ID_EQ, COLOR_ID_EQ + COLOR_ID_EQ_LEN);
    buf.insert(buf.end(), color_set_id, color_set_id + color_set_id_len);
    buf.push_back('\n');

    buf.insert(buf.end(), unitig_label, unitig_label + unitig_label_length);
    buf.push_back('\n');

    unitigs_out.write(buf.data(), buf.size());
}

}