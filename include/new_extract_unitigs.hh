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

void write_unitig(
    const char* unitig_id_chars, int64_t unitig_id_length, 
    const char* unitig_label, int64_t unitig_label_length, 
    const char* color_set_id, int64_t color_set_id_len, 
    ParallelOutputWriter& unitigs_out);

void write_colorset(int64_t color_set_id, const vector<int64_t>& colors, ParallelOutputWriter& out);

// Splits unitigs by colorset runs
// Returns the DBG nodes that were visited
// This function needs to be thread-safe
template<typename coloring_t>
vector<DBG::Node> process_unitig_from(const DBG& dbg, const coloring_t& coloring, DBG::Node v, ParallelOutputWriter& unitigs_out, ParallelOutputWriter& colors_out) {

    vector<DBG::Node> nodes;
    vector<int64_t> subunitig_ends; // Unitigs broken by color set runs. Exclusive endpoints
    vector<char> label;
    std::tie(nodes, label) = walk_unitig_from(dbg, v);
    check_true(nodes.size() > 0, "BUG: empty unitig");

    int64_t color_set_id = coloring.get_color_set_id(nodes.back().id);
    vector<int64_t> color_set_ids;

    for(int64_t pos = (int64_t)nodes.size()-1; pos >= 0; pos--){
        const DBG::Node& u = nodes[pos];

        // Color sets can change only at core k-mers. Otherwise the color set is the same as that
        // of the successor in the DBG.
        int64_t new_color_set_id = coloring.is_core_kmer(u.id) ? coloring.get_color_set_id(u.id) : color_set_id;
        if(pos == (int64_t)nodes.size()-1 || new_color_set_id != color_set_id) {
            subunitig_ends.push_back(pos+1); // Start a new subunitig
            color_set_ids.push_back(new_color_set_id);
        }
        color_set_id = new_color_set_id;
    }

    subunitig_ends.push_back(0); // Sentinel
    color_set_ids.push_back(-1); // Sentinel
    std::reverse(subunitig_ends.begin(), subunitig_ends.end());
    std::reverse(color_set_ids.begin(), color_set_ids.end());

    char unitig_id_buf[32]; // Enough space to encode 64-bit integers in ascii
    char color_set_id_buf[32]; // Enough space to encode 64-bit integers in ascii
    for(int64_t i = 1; i < subunitig_ends.size(); i++) {

        int64_t unitig_id = nodes[subunitig_ends[i-1]].id; // Unitig id is the colex rank of the first k-mer of the subunitig
        int64_t unitig_id_string_len = fast_int_to_string(unitig_id, unitig_id_buf);

        int64_t color_set_id = color_set_ids[i]; 
        int64_t color_set_id_string_len = fast_int_to_string(color_set_id, color_set_id_buf);

        int64_t len = subunitig_ends[i] - subunitig_ends[i-1]; // Length in nodes
        int64_t string_len = len + (dbg.get_k() - 1); // Length of the string label

        write_unitig(unitig_id_buf, unitig_id_string_len, label.data() + subunitig_ends[i-1], string_len, color_set_id_buf, color_set_id_string_len, unitigs_out);

    }

    return nodes;

}

template<typename coloring_t>
void write_distinct_color_sets(const coloring_t& coloring, ParallelOutputWriter& out){

    int64_t n_sets = coloring.number_of_distinct_color_sets();

    vector<int64_t> color_set_buf;
    for(int64_t color_set_id = 0; color_set_id < n_sets; color_set_id++){
        typename coloring_t::colorset_view_type color_set = coloring.get_color_set_by_color_set_id(color_set_id);
        color_set.push_colors_to_vector(color_set_buf);
        
        write_colorset(color_set_id, color_set_buf, out);
        color_set_buf.clear();
    }
}

} // End of namespace

template<typename coloring_t>
void new_extract_unitigs(int64_t n_threads, const DBG& dbg, string unitigs_outfile, optional<coloring_t*> coloring,
                         optional<string> colorsets_outfile){}; // TODO: remove

/*

salmonella_10.metadata.txt:

num_references=10
num_unitigs=86630
num_color_classes=171

salmonella_10.unitigs.fa:

> unitig_id=0 color_id=0
GATTGAGCACCAACTGCGAGAATCAGGTGTTGAAGAGCAAGGGCGTGTGTTTATCGAAAAAGCTATTGAGCAGCCGCTTGATCCACAA
> unitig_id=1 color_id=0
GAAATTTAACGGCTGTTTTTCCGGCCAGATGTTATGTCTGGCTGGTTTTATTGTTTTGATTTTAAAGGAATTTACAGTGAATAAATGGCGTAACCCCACTGGGTGGTTATGTGCGGTAGCTATGCCTTTTG
> unitig_id=2 color_id=0
GCGCTGAACATCAGCGCCTTTCTGCGACAGCTCAATCATGCATTCGCCAATCACGGCAATC
...

salmonella_10.colors.txt:

color_id=0 size=4 0 3 7 8
color_id=1 size=1 8
color_id=2 size=10 0 1 2 3 4 5 6 7 8 9
color_id=3 size=6 1 2 4 5 6 9
...

*/

template<typename coloring_t>
void dump_index(int64_t n_threads, const DBG& dbg, coloring_t& coloring, optional<string> unitigs_outfile, optional<string> colorsets_outfile, optional<string> metadata_outfile) {

    using namespace new_extract_unitigs_internals;

    // Output writers. Won't write anything if the given optional filename is null.
    ParallelOutputWriter unitigs_out(unitigs_outfile); // Thread-safe output writer
    ParallelOutputWriter colors_out(colorsets_outfile); // Thread-safe output writer
    ParallelOutputWriter metadata_out(metadata_outfile); // Thread-safe output writer

    if(unitigs_out.enabled) {
        vector<bool> visited(dbg.number_of_sets_in_sbwt());

        write_log("Constructing acyclic unitigs", LogLevel::MAJOR);
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
        write_log("Constructing cyclic unitigs", LogLevel::MAJOR);
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
    }

    if(colors_out.enabled){
        write_log("Writing color sets", LogLevel::MAJOR);
        write_distinct_color_sets(coloring, colors_out);
        colors_out.flush();
    }

    metadata_out.flush();

    write_log("Done writing unitigs", LogLevel::MAJOR);

}