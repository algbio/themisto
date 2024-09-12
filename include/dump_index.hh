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

// Break unitigs by by color set runs. Exclusive endpoints.
template<typename coloring_t>
pair<vector<int64_t>, vector<int64_t>> get_subunitig_ends_and_color_set_ids(const coloring_t& coloring, const vector<DBG::Node>& nodes) {
    int64_t color_set_id = coloring.get_color_set_id(nodes.back().id);
    vector<int64_t> color_set_ids;

    vector<int64_t> subunitig_ends;

    vector<int64_t> prev_color_set;
    coloring.get_color_set_by_color_set_id(color_set_id).push_colors_to_vector(prev_color_set);

    vector<int64_t> cur_color_set;

    for(int64_t pos = (int64_t)nodes.size()-1; pos >= 0; pos--){
        const DBG::Node& u = nodes[pos];

        // Color sets can change only at core k-mers. Otherwise the color set is the same as that
        // of the successor in the DBG.
        int64_t new_color_set_id = coloring.is_core_kmer(u.id) ? coloring.get_color_set_id(u.id) : color_set_id;

        bool color_set_changed = false;

        if((pos == (int64_t)nodes.size()-1)){ // Initialize the color set for the first one
            coloring.get_color_set_by_color_set_id(new_color_set_id).push_colors_to_vector(cur_color_set);
            color_set_changed = true; // The first color set considered a new one
        }

        if(new_color_set_id != color_set_id) {
            // Color set may have changed. It can sometimes happen that two color sets with different ids
            // are actually the same. This is not ideal but GGCAT seems to sometimes do this during unitig
            // construction. We could add a step to eliminate duplicate color sets, in which case it would
            // be enough to compare ids, but for now, we compare the actual color sets.
            cur_color_set.clear();
            coloring.get_color_set_by_color_set_id(new_color_set_id).push_colors_to_vector(cur_color_set);
            color_set_changed |= prev_color_set != cur_color_set;
        }

        if(color_set_changed) {
            subunitig_ends.push_back(pos+1); // Start a new subunitig
            color_set_ids.push_back(new_color_set_id);

            prev_color_set.clear();
            std::swap(prev_color_set, cur_color_set);
        }

        color_set_id = new_color_set_id;
        
    }

    subunitig_ends.push_back(0); // Sentinel
    color_set_ids.push_back(-1); // Sentinel
    std::reverse(subunitig_ends.begin(), subunitig_ends.end());
    std::reverse(color_set_ids.begin(), color_set_ids.end());

    return {subunitig_ends, color_set_ids};

}

struct UnitigBothWays {
    vector<DBG::Node> fw_nodes;
    string fw_label;
    vector<DBG::Node> rc_nodes;
    string rc_label;

    pair<vector<DBG::Node>&, string> canonical() {
        if(fw_label < rc_label) return {fw_nodes, fw_label};
        else return {rc_nodes, rc_label};
    }

};

UnitigBothWays get_both_unitig_orientations(const DBG& dbg, DBG::Node start_node);

// Splits unitigs by colorset runs.
// This function needs to be thread-safe
template<typename coloring_t>
vector<UnitigBothWays> process_unitig_from(const DBG& dbg, const coloring_t& coloring, DBG::Node v, ParallelOutputWriter& unitigs_out) {

    UnitigBothWays unitig_both_ways = get_both_unitig_orientations(dbg, v);

    //const vector<DBG::Node>& nodes = unitig_both_ways.canonical().first; 
    //const string& label = unitig_both_ways.canonical().second; 

    vector<int64_t> subunitig_ends, color_set_ids;
    tie(subunitig_ends, color_set_ids) = get_subunitig_ends_and_color_set_ids(coloring, unitig_both_ways.fw_nodes);

    vector<UnitigBothWays> subunitigs;
    //char unitig_id_buf[32]; // Enough space to encode 64-bit integers in ascii
    //char color_set_id_buf[32]; // Enough space to encode 64-bit integers in ascii
    for(int64_t i = 1; i < subunitig_ends.size(); i++) {

        /*
        int64_t unitig_id = nodes[subunitig_ends[i-1]].id; // Unitig id is the colex rank of the first k-mer of the subunitig. This is before canonicalizatoin but that's ok.
        int64_t unitig_id_string_len = fast_int_to_string(unitig_id, unitig_id_buf);

        int64_t color_set_id = color_set_ids[i]; 
        int64_t color_set_id_string_len = fast_int_to_string(color_set_id, color_set_id_buf);
        */

        int64_t node_len = subunitig_ends[i] - subunitig_ends[i-1]; // Length in nodes
        int64_t string_len = node_len + (dbg.get_k() - 1); // Length of the string label

        int64_t start = subunitig_ends[i-1]; // Start node index

        string fw_label = unitig_both_ways.fw_label.substr(start, string_len);
        vector<int64_t> fw_nodes(unitig_both_ways.fw_nodes.begin() + start, unitig_both_ways.fw_nodes.begin() + start + node_len); 

        int64_t rc_start = unitig_both_ways.fw_nodes.size() - 1 - (start + node_len - 1);

        string rc_label = unitig_both_ways.rc_label.substr(rc_start, string_len);
        vector<int64_t> rc_nodes(unitig_both_ways.rc_nodes.begin() + rc_start, unitig_both_ways.rc_nodes.begin() + rc_start + node_len); 

        UnitigBothWays subunitig = {fw_nodes, fw_label, rc_nodes, rc_label};
        subunitigs.push_back(subunitig);

        //write_unitig(unitig_id_buf, unitig_id_string_len, label.data() + subunitig_ends[i-1], string_len, color_set_id_buf, color_set_id_string_len, unitigs_out);
    }

    return subunitigs;

}

template<typename coloring_t>
void write_distinct_color_sets(const coloring_t& coloring, ParallelOutputWriter& out){

    int64_t n_sets = coloring.number_of_distinct_color_sets();

    vector<int64_t> color_set_buf;
    Progress_printer pp(n_sets, 100);
    for(int64_t color_set_id = 0; color_set_id < n_sets; color_set_id++){
        typename coloring_t::colorset_view_type color_set = coloring.get_color_set_by_color_set_id(color_set_id);
        color_set.push_colors_to_vector(color_set_buf);
        
        write_colorset(color_set_id, color_set_buf, out);
        color_set_buf.clear();

        pp.job_done();
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
void dump_index(int64_t n_threads, const DBG& dbg, coloring_t& coloring, optional<string> unitigs_outfile, optional<string> colorsets_outfile, optional<string> metadata_outfile, optional<string> sbwt_outfile) {

    using namespace new_extract_unitigs_internals;

    // Output writers. Won't write anything if the given optional filename is null.
    ParallelOutputWriter unitigs_out(unitigs_outfile); // Thread-safe output writer
    ParallelOutputWriter colors_out(colorsets_outfile); // Thread-safe output writer
    ParallelOutputWriter metadata_out(metadata_outfile); // Thread-safe output writer

    if(unitigs_out.enabled) {
        vector<bool> visited(dbg.number_of_sets_in_sbwt());

        int64_t num_unitigs = 0;
        write_log("Constructing acyclic unitigs", LogLevel::MAJOR);
        Progress_printer pp_acyclic(dbg.number_of_kmers(), 100);
        #pragma omp parallel for num_threads (n_threads)
        for(int64_t colex = 0; colex < dbg.number_of_sets_in_sbwt(); colex++){
            #pragma omp critical // TODO: this might cause a lot of parallel contention
            {
                pp_acyclic.job_done();
            }

            if(dbg.is_dummy_colex_position(colex)) continue;

            DBG::Node v(colex);
            if(!is_first_kmer_of_unitig(dbg, v)) continue;

            UnitigBothWays unitig = process_unitig_from(dbg, coloring, v, unitigs_out);

            #pragma omp critical // Modifying shared data
            {
                for(DBG::Node u : unitig.fw_nodes){
                    assert(!visited[u.id]);
                    visited[u.id] = true;
                }
                for(DBG::Node u : unitig.rc_nodes){
                    assert(!visited[u.id]);
                    visited[u.id] = true;
                }
                num_unitigs++;
            }
        }

        // Only disjoint cyclic unitigs remain
        write_log("Constructing cyclic unitigs", LogLevel::MAJOR);
        Progress_printer pp_cyclic(dbg.number_of_kmers(), 100);
        for(DBG::Node v : dbg.all_nodes()) {
            pp_cyclic.job_done();
            if(visited[v.id]) continue;
            UnitigBothWays unitig = process_unitig_from(dbg, coloring, v, unitigs_out);
            for(DBG::Node u : unitig.fw_nodes){
                assert(!visited[u.id]);
                visited[u.id] = true;
            }
            for(DBG::Node u : unitig.rc_nodes){
                assert(!visited[u.id]);
                visited[u.id] = true;
            }
            num_unitigs++;
        }

        unitigs_out.flush();

        metadata_out.write("num_unitigs=");
        metadata_out.write(to_string(num_unitigs));
        metadata_out.write("\n");
    }

    if(colors_out.enabled){
        write_log("Writing color sets", LogLevel::MAJOR);
        write_distinct_color_sets(coloring, colors_out);
        colors_out.flush();

        metadata_out.write("num_colors=");
        metadata_out.write(to_string(coloring.largest_color()+1));
        metadata_out.write("\n");

        metadata_out.write("num_color_sets=");
        metadata_out.write(to_string(coloring.number_of_distinct_color_sets()));
        metadata_out.write("\n");

        metadata_out.write("k=");
        metadata_out.write(to_string(dbg.get_k()));
        metadata_out.write("\n");
    }

    if(sbwt_outfile.has_value()){
        write_log("Writing SBWT", LogLevel::MAJOR);
        const plain_matrix_sbwt_t& sbwt = coloring.get_sbwt();
        throwing_ofstream out(sbwt_outfile.value());
        sbwt.ascii_export_metadata(out.stream);
        sbwt.ascii_export_sets(out.stream);
    }

    metadata_out.flush();

    write_log("Index dump finished", LogLevel::MAJOR);

}
