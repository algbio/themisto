#pragma once

#include <iostream>
#include "DBG.hh"
#include "globals.hh"
#include "coloring/Coloring.hh"

bool is_first_kmer_of_unitig(const DBG& dbg, const DBG::Node& node) {
    int64_t indeg = dbg.indegree(node);
    if(indeg != 1){
        return true;
    } else {
        DBG::Node u = dbg.pred(node);
        return dbg.outdegree(u) > 1;
    }
}

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

template<typename coloring_t>
void new_extract_unitigs(const DBG& dbg, const coloring_t& coloring, ostream& unitigs_out,
                         bool split_by_colorset_runs, ostream& colorsets_out,
                         ostream& gfa_out, int64_t min_colors = 0) {

    int64_t unitig_id = 0;
    vector<bool> visited(dbg.SBWT->n_nodes());

    write_log("Listing acyclic unitigs", LogLevel::MAJOR);
    for(DBG::Node v : dbg.all_nodes()){
        if(!is_first_kmer_of_unitig(dbg, v)) continue;
        vector<int64_t> nodes;
        vector<char> label;
        std::tie(nodes,label) = walk_unitig_from(dbg, v);
        for(int64_t u : nodes){
            assert(!visited[u]);
            visited[u] = true;
        }
        label.push(0); // Make it a C string
        unitigs_out << ">" << unitig_id++ << label.data() << "\n";

        // TODO: break by color sets
    }

    // Only disjoint cyclic unitigs remain
    write_log("Listing cyclic unitigs", LogLevel::MAJOR);

    for(DBG::Node v : dbg.all_nodes()) {
        if(visited[v]) continue;

        vector<int64_t> nodes;
        vector<char> label;
        std::tie(nodes,label) = walk_unitig_from(dbg, v);

        for(int64_t u : nodes){
            assert(!visited[u]);
            visited[u] = true;
        }

        label.push(0); // Make it a C string
        unitigs_out << ">" << unitig_id++ << label.data() << "\n";

        // TODO: break by color sets
    }

    write_log("Done writing unitigs", LogLevel::MAJOR);

}