#pragma once

#include <iostream>
#include "DBG.hh"
#include "globals.hh"
#include "coloring/Coloring.hh"

bool is_first_kmer_of_unitig(const DBG& dbg, const DBG::Node& node); 

// Returns the sequence of nodes and the label of the unitig
pair<vector<DBG::Node>, vector<char>> walk_unitig_from(const DBG& dbg, DBG::Node v);

// If coloring is given, splits unitigs by colorset runs
template<typename coloring_t, typename out_stream_t>
void process_unitig_from(const DBG& dbg, optional<coloring_t*> coloring, DBG::Node v, vector<bool>& visited, out_stream_t& unitigs_out, int64_t unitig_id) {

    vector<DBG::Node> nodes;
    vector<int64_t> subunitig_ends; // Unitigs broken by color set runs. Exclusive endpoints
    vector<char> label;
    std::tie(nodes, label) = walk_unitig_from(dbg, v);
    check_true(nodes.size() > 0, "BUG: empty unitig");

    int64_t color_set_id = 0;
    if(coloring.has_value()) color_set_id = coloring.value()->get_color_set_id(nodes.back().id);

    for(int64_t pos = (int64_t)nodes.size()-1; pos >= 0; pos--){
        const DBG::Node& u = nodes[pos];

        if(coloring.has_value()){
            // Color sets can change only at core k-mers. Otherwise the color set is the same as that
            // of the successor in the DBG.
            int64_t new_color_set_id = coloring.value()->is_core_kmer(u.id) ? coloring.value()->get_color_set_id(u.id) : color_set_id;
            if(pos == (int64_t)nodes.size()-1 || new_color_set_id != color_set_id) {
                subunitig_ends.push_back(pos+1); // Start a new subunitig
            }
            color_set_id = new_color_set_id;
        }

        assert(!visited[u.id]);
        visited[u.id] = true;
    }

    if(coloring.has_value()){
        subunitig_ends.push_back(0); // Sentinel
        std::reverse(subunitig_ends.begin(), subunitig_ends.end());
    } else {
        subunitig_ends = {0, nodes.size()}; // One big unitig
    }

    char unitig_id_buf[32]; // Enough space to encode 64-bit integers in ascii
    for(int64_t i = 1; i < subunitig_ends.size(); i++) {
        int64_t unitig_id_string_len = fast_int_to_string(unitig_id++, unitig_id_buf);
        unitigs_out.write(">", 1);
        unitigs_out.write(unitig_id_buf, unitig_id_string_len);
        unitigs_out.write("\n", 1);

        int64_t len = subunitig_ends[i] - subunitig_ends[i-1]; // Length in nodes
        int64_t string_len = len + (dbg.get_k() - 1); // Length of the string label
        unitigs_out.write(label.data() + subunitig_ends[i-1], string_len);
        unitigs_out.write("\n", 1);
    }

}

template<typename coloring_t, typename out_stream_t>
void new_extract_unitigs(const DBG& dbg, out_stream_t& unitigs_out, optional<coloring_t*> coloring,
                         optional<out_stream_t*> colorsets_out, 
                         int64_t min_colors = 0) {

    // TODO: min_colors
    // TODO: GFA

    int64_t unitig_id = 0;
    vector<bool> visited(dbg.number_of_sets_in_sbwt());

    write_log("Listing acyclic unitigs", LogLevel::MAJOR);
    Progress_printer pp_acyclic(dbg.number_of_kmers(), 100);
    for(DBG::Node v : dbg.all_nodes()){
        pp_acyclic.job_done();
        if(!is_first_kmer_of_unitig(dbg, v)) continue;
        process_unitig_from(dbg, coloring, v, visited, unitigs_out, unitig_id++);
    }

    // Only disjoint cyclic unitigs remain
    write_log("Listing cyclic unitigs", LogLevel::MAJOR);
    Progress_printer pp_cyclic(dbg.number_of_kmers(), 100);
    for(DBG::Node v : dbg.all_nodes()) {
        pp_cyclic.job_done();
        if(visited[v.id]) continue;
        process_unitig_from(dbg, coloring, v, visited, unitigs_out, unitig_id++);
    }

    unitigs_out.flush();
    write_log("Done writing unitigs", LogLevel::MAJOR);

}