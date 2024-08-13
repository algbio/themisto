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

static const char* COLOR_SET_ID_EQ = "color_set_id=";
static const int64_t COLOR_SET_ID_EQ_LEN = 13;

static const char* SIZE_EQ = "size=";
static const int64_t SIZE_EQ_LEN= 5;

void write_unitig(
    const char* unitig_id_chars, int64_t unitig_id_length, 
    const char* unitig_label, int64_t unitig_label_length, 
    const char* color_set_id, int64_t color_set_id_len, 
    ParallelOutputWriter& unitigs_out){

    vector<char> buf; // ASCII data. Needs to be written in one go to avoid interleaved writing with other threads!

    buf.push_back('>');
    buf.push_back(' ');
    buf.insert(buf.end(), UNITIG_ID_EQ, UNITIG_ID_EQ + UNITIG_ID_EQ_LEN);
    buf.insert(buf.end(), unitig_id_chars, unitig_id_chars + unitig_id_length); // Push the unitig id
    buf.push_back(' ');
    buf.insert(buf.end(), COLOR_SET_ID_EQ, COLOR_SET_ID_EQ + COLOR_SET_ID_EQ_LEN); // Push color set id
    buf.insert(buf.end(), color_set_id, color_set_id + color_set_id_len);
    buf.push_back('\n');

    buf.insert(buf.end(), unitig_label, unitig_label + unitig_label_length);
    buf.push_back('\n');

    unitigs_out.write(buf.data(), buf.size());
}

void write_colorset(int64_t color_set_id, const vector<int64_t>& colors, ParallelOutputWriter& out){

    // Write a line like:
    // color_id=0 size=4 0 3 7 8

    vector<char> buf; // ASCII line that will be written. Needs to be written in one go to avoid interleaved writing with other threads!
    char int_buf[32]; // Enough space to encode 64-bit integers in ascii
    int64_t int_buf_len = 0;

    // Write "color_set_id="
    buf.insert(buf.end(), COLOR_SET_ID_EQ, COLOR_SET_ID_EQ + COLOR_SET_ID_EQ_LEN);

    // Write color id integer
    int_buf_len = fast_int_to_string(color_set_id, int_buf); 
    buf.insert(buf.end(), int_buf, int_buf + int_buf_len);

    buf.push_back(' ');

    // Write "size="
    buf.insert(buf.end(), SIZE_EQ, SIZE_EQ + SIZE_EQ_LEN);

    // Write the size integer
    int_buf_len = fast_int_to_string(colors.size(), int_buf); 
    buf.insert(buf.end(), int_buf, int_buf + int_buf_len);

    // Write color ids
    for(int64_t color : colors){
        buf.push_back(' ');
        int_buf_len = fast_int_to_string(color, int_buf); 
        buf.insert(buf.end(), int_buf, int_buf + int_buf_len);
    }

    buf.push_back('\n');

    out.write(buf.data(), buf.size());
}


}