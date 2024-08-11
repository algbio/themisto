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
void process_unitig_from(const DBG& dbg, const coloring_t& coloring, DBG::Node v, vector<bool>& visited, ostream& unitigs_out, int64_t unitig_id, bool split_by_colorset_runs) {
    vector<DBG::Node> nodes;
    vector<int64_t> subunitig_ends; // Unitigs broken by color set runs. Exclusive endpoints
    vector<char> label;
    std::tie(nodes, label) = walk_unitig_from(dbg, v);
    check_true(nodes.size() > 0, "BUG: empty unitig");

    int64_t color_set_id = coloring.get_color_set_id(nodes.back().id);

    for(int64_t pos = (int64_t)nodes.size()-1; pos >= 0; pos--){
        const DBG::Node& u = nodes[pos];

        if(split_by_colorset_runs){
            // Color sets can change only at core k-mers. Otherwise the color set is the same as that
            // of the successor in the DBG.
            int64_t new_color_set_id = coloring.is_core_kmer(u.id) ? coloring.get_color_set_id(u.id) : color_set_id;
            if(pos == (int64_t)nodes.size()-1 || new_color_set_id != color_set_id) {
                subunitig_ends.push_back(pos+1); // Start a new subunitig
            }
            color_set_id = new_color_set_id;
        }

        assert(!visited[u.id]);
        visited[u.id] = true;
    }

    if(split_by_colorset_runs){
        subunitig_ends.push_back(0); // Sentinel
        std::reverse(subunitig_ends.begin(), subunitig_ends.end());
    } else {
        subunitig_ends = {0, nodes.size()}; // One big unitig
    }
    for(int64_t i = 1; i < subunitig_ends.size(); i++) {
        int64_t len = subunitig_ends[i] - subunitig_ends[i-1]; // Length in nodes
        int64_t string_len = len + (dbg.get_k() - 1); // Length of the string label
        unitigs_out << ">" << unitig_id << "\n";
        unitigs_out.write(label.data() + subunitig_ends[i-1], string_len);
        unitigs_out.write("\n", 1);
    }

}

template<typename coloring_t>
void new_extract_unitigs(const DBG& dbg, const coloring_t& coloring, ostream& unitigs_out,
                         bool split_by_colorset_runs, ostream& colorsets_out,
                         ostream& gfa_out, int64_t min_colors = 0) {

    int64_t unitig_id = 0;
    vector<bool> visited(dbg.number_of_sets_in_sbwt());

    write_log("Listing acyclic unitigs", LogLevel::MAJOR);
    for(DBG::Node v : dbg.all_nodes()){
        if(!is_first_kmer_of_unitig(dbg, v)) continue;
        process_unitig_from(dbg, coloring, v, visited, unitigs_out, unitig_id++, split_by_colorset_runs);
    }

    // Only disjoint cyclic unitigs remain
    write_log("Listing cyclic unitigs", LogLevel::MAJOR);

    for(DBG::Node v : dbg.all_nodes()) {
        if(visited[v.id]) continue;
        process_unitig_from(dbg, coloring, v, visited, unitigs_out, unitig_id++, split_by_colorset_runs);
    }

    unitigs_out.flush();
    write_log("Done writing unitigs", LogLevel::MAJOR);

}