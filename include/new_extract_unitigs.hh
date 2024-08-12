#pragma once

#include <iostream>
#include "DBG.hh"
#include "globals.hh"
#include "coloring/Coloring.hh"
#include "WorkDispatcher.hh"

namespace new_extract_unitigs_internals { // Wrap in namespace to avoid name collisions

bool is_first_kmer_of_unitig(const DBG& dbg, const DBG::Node& node); 

// Returns the sequence of nodes and the label of the unitig
pair<vector<DBG::Node>, vector<char>> walk_unitig_from(const DBG& dbg, DBG::Node v);

void write_unitig(const char* unitig_id_chars, int64_t unitig_id_length, const char* unitig_label, int64_t unitig_label_length, ParallelOutputWriter& unitigs_out);

template<typename colorset_view_t> 
void write_colorset(const char* unitig_id_chars, int64_t unitig_id_length, colorset_view_t& colorset, ParallelOutputWriter& colors_out){
    vector<char> color_set_line; // ASCII line that will be written. Needs to be written in one go to avoid interleaved writing with other threads!
    char color_id_buf[32]; // Enough space to encode 64-bit integers in ascii

    color_set_line.insert(color_set_line.end(), unitig_id_chars, unitig_id_chars + unitig_id_length); // Push the unitig id
    for (int64_t color : colorset.get_colors_as_vector()) {
        int64_t color_id_string_len = fast_int_to_string(color, color_id_buf);
        color_set_line.push_back(' ');
        color_set_line.insert(color_set_line.end(), color_id_buf, color_id_buf + color_id_string_len);
    }
    color_set_line.push_back('\n');
    colors_out.write(color_set_line.data(), color_set_line.size());
}


// If coloring is given, splits unitigs by colorset runs
// Returns the DBG nodes that were visited
// This function needs to be thread-safe
template<typename coloring_t>
vector<DBG::Node> process_unitig_from(const DBG& dbg, const optional<coloring_t*> coloring, DBG::Node v, ParallelOutputWriter& unitigs_out, optional<ParallelOutputWriter>& colors_out) {

    vector<DBG::Node> nodes;
    vector<int64_t> subunitig_ends; // Unitigs broken by color set runs. Exclusive endpoints
    vector<char> label;
    std::tie(nodes, label) = walk_unitig_from(dbg, v);
    check_true(nodes.size() > 0, "BUG: empty unitig");

    int64_t color_set_id = 0;
    if(coloring.has_value()) color_set_id = coloring.value()->get_color_set_id(nodes.back().id);
    vector<int64_t> color_set_ids;

    for(int64_t pos = (int64_t)nodes.size()-1; pos >= 0; pos--){
        const DBG::Node& u = nodes[pos];

        if(coloring.has_value()){
            // Color sets can change only at core k-mers. Otherwise the color set is the same as that
            // of the successor in the DBG.
            int64_t new_color_set_id = coloring.value()->is_core_kmer(u.id) ? coloring.value()->get_color_set_id(u.id) : color_set_id;
            if(pos == (int64_t)nodes.size()-1 || new_color_set_id != color_set_id) {
                subunitig_ends.push_back(pos+1); // Start a new subunitig
                color_set_ids.push_back(new_color_set_id);
            }
            color_set_id = new_color_set_id;
        }

    }

    if(coloring.has_value()){
        subunitig_ends.push_back(0); // Sentinel
        color_set_ids.push_back(-1); // Sentinel
        std::reverse(subunitig_ends.begin(), subunitig_ends.end());
        std::reverse(color_set_ids.begin(), color_set_ids.end());
    } else {
        subunitig_ends = {0, nodes.size()}; // One big unitig
    }

    char unitig_id_buf[32]; // Enough space to encode 64-bit integers in ascii
    for(int64_t i = 1; i < subunitig_ends.size(); i++) {

        int64_t unitig_id = nodes[subunitig_ends[i-1]].id; // Unitig id is the colex rank of the first k-mer of the subunitig
        int64_t unitig_id_string_len = fast_int_to_string(unitig_id, unitig_id_buf);
        int64_t len = subunitig_ends[i] - subunitig_ends[i-1]; // Length in nodes
        int64_t string_len = len + (dbg.get_k() - 1); // Length of the string label

        write_unitig(unitig_id_buf, unitig_id_string_len, label.data() + subunitig_ends[i-1], string_len, unitigs_out);

        if(colors_out.has_value()) { // Write color set
            typename coloring_t::colorset_view_type colorset = (*coloring)->get_color_set_by_color_set_id(color_set_ids[i]);
            write_colorset(unitig_id_buf, unitig_id_string_len, colorset, *colors_out);
        }
    }

    return nodes;

}

} // End of namespace

template<typename coloring_t>
void new_extract_unitigs(int64_t n_threads, const DBG& dbg, string unitigs_outfile, optional<coloring_t*> coloring,
                         optional<string> colorsets_outfile,
                         int64_t min_colors = 0) {

    using namespace new_extract_unitigs_internals;

    // TODO: min_colors

    ParallelOutputWriter unitigs_out(unitigs_outfile); // Thread-safe output writer

    optional<ParallelOutputWriter> colors_out;
    if(colorsets_outfile.has_value())
        colors_out.emplace(colorsets_outfile.value());

    vector<bool> visited(dbg.number_of_sets_in_sbwt());

    write_log("Listing acyclic unitigs", LogLevel::MAJOR);
    Progress_printer pp_acyclic(dbg.number_of_kmers(), 100);
    #pragma omp parallel for num_threads (n_threads)
    for(int64_t colex = 0; colex < dbg.number_of_sets_in_sbwt(); colex++){
        #pragma omp critical
        {
            pp_acyclic.job_done();
        }

        if(dbg.is_dummy_colex_position(colex)) continue;

        DBG::Node v(colex);
        if(!is_first_kmer_of_unitig(dbg, v)) continue;

        vector<DBG::Node> nodes = process_unitig_from(dbg, coloring, v, unitigs_out, colors_out);

        #pragma omp critical // Modifying shared data
        {
            for(DBG::Node u : nodes){
                assert(!visited[u.id]);
                visited[u.id] = true;
            }
        }
    }

    // Only disjoint cyclic unitigs remain
    write_log("Listing cyclic unitigs", LogLevel::MAJOR);
    Progress_printer pp_cyclic(dbg.number_of_kmers(), 100);
    for(DBG::Node v : dbg.all_nodes()) {
        pp_cyclic.job_done();
        if(visited[v.id]) continue;
        vector<DBG::Node> nodes = process_unitig_from(dbg, coloring, v, unitigs_out, colors_out);
        for(DBG::Node u : nodes){
            assert(!visited[u.id]);
            visited[u.id] = true;
        }
    }

    unitigs_out.flush();
    if(colors_out.has_value()) colors_out->flush();

    write_log("Done writing unitigs", LogLevel::MAJOR);

}