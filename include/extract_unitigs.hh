#pragma once

#include <iostream>

#include "DBG.hh"
#include "globals.hh"
#include "new_coloring.hh"

using namespace std;

class UnitigExtractor {
   private:
    struct Unitig {
        vector<DBG::Node> nodes;
        vector<int64_t>
            links;  // Outgoing edges from this unitig, as a list of unitig ids
        int64_t id;
    };

    struct Colored_Unitig {
        vector<DBG::Node> nodes;
        vector<uint32_t> colorset;
        vector<int64_t>
            links;  // Outgoing edges from this unitig, as a list of unitig ids
        int64_t id;
    };

    // Assumes that dummy nodes are always visited
    Unitig get_unitig_containing_node(DBG::Node v, const DBG& dbg, unordered_map<DBG::Node, bool>& visited) {
        // TODO: visited hash map is slow.
        vector<DBG::Node> backward, forward;

        visited[v] = true;

        // Walk backward until (indegree >= 2), or (predecessor is visited or
        // has outdegree >= 2)
        DBG::Node u = v;
        while (true) {
            int64_t indeg = dbg.indegree(u);
            if (indeg >= 2 || indeg == 0) break;
            DBG::Node pred = dbg.pred(u);
            if (pred.id == -1) break;
            if (visited[pred]) break;
            if (dbg.outdegree(pred) >= 2) break;
            u = pred;
            visited[u] = true;
            backward.push_back(u);
        }

        // Walk forward until (outdegree >= 2) or (successor is visited or has
        // indegree >= 2)
        u = v;
        while (true) {
            int64_t outdeg = dbg.outdegree(u);
            if (outdeg >= 2 || outdeg == 0) break;
            DBG::Node succ = dbg.succ(u);
            if (succ.id == -1) break;
            if (visited[succ]) break;
            if (dbg.outdegree(succ) >= 2) break;
            u = succ;
            visited[u] = true;
            forward.push_back(u);
        }

        // Collect the unitig path
        Unitig U;
        for (int64_t i = (int64_t)backward.size() - 1; i >= 0; i--)
            U.nodes.push_back(backward[i]);
        U.nodes.push_back(v);
        for (DBG::Node w : forward) U.nodes.push_back(w);

        // Compute links
            
        for(DBG::Edge edge : dbg.outedges(U.nodes.back())){
            DBG::Node destination = edge.dest;
            U.links.push_back(destination.id);
        }

        U.id = U.nodes[0].id;
        return U;
    }

    vector<Colored_Unitig> split_to_colorset_runs(Unitig& U, const Coloring& coloring) {
        vector<Colored_Unitig> colored_unitigs;
        auto get_id = [&](int64_t node) {
            return coloring.get_color_set_id(node);
        };
        for (int64_t run_start = 0; run_start < U.nodes.size(); run_start++) {
            int64_t run_end = run_start;
            while (run_end < U.nodes.size() - 1 &&
                   get_id(U.nodes[run_start].id) == get_id(U.nodes[run_end + 1].id)) {
                run_end++;
            }

            Colored_Unitig colored;
            for (int64_t i = run_start; i <= run_end; i++)
                colored.nodes.push_back(U.nodes[i]);
            colored.colorset = coloring.get_color_set_as_vector(U.nodes[run_start].id);
            colored.id = U.nodes[run_start].id;

            colored_unitigs.push_back(colored);  // Links are added later

            run_start = run_end;
        }

        // Linkage
        for (int64_t i = 0; i < colored_unitigs.size(); i++) {
            if (i < (int64_t)colored_unitigs.size() - 1)  // Not last
                colored_unitigs[i].links.push_back(
                    colored_unitigs[i + 1].id);  // Chain unitigs
            else                                 // Last
                colored_unitigs[i].links =
                    U.links;  // The outgoing links of the original unitig
        }

        return colored_unitigs;
    }

    // Extracts maximal unitigs that have at least n_min_colors colors
    vector<Colored_Unitig> split_by_colorset_size(const int64_t n_min_colors,
                                                  Unitig& U,
                                                  const Coloring& coloring) {
        vector<Colored_Unitig> colored_unitigs;
        auto get_n_colors = [&](DBG::Node node) {
            return coloring.get_color_set(node.id).size();
        };

        int64_t run_start = 0;
        while (run_start < U.nodes.size()) {
            int64_t n_colors_for_node = get_n_colors(U.nodes[run_start]);

            if (n_colors_for_node >= n_min_colors) {
                int64_t run_end = run_start;
                while (run_end < U.nodes.size() - 1 &&
                       get_n_colors(U.nodes[run_end + 1]) >= n_min_colors) {
                    run_end++;
                }

                Colored_Unitig colored;
                for (int64_t i = run_start; i <= run_end; i++)
                    colored.nodes.push_back(U.nodes[i]);

                colored.id = U.nodes[run_start].id;

                colored_unitigs.push_back(colored);  // Links are added later

                run_start = run_end;
            }
            ++run_start;
        }

        // Linkage
        for (int64_t i = 0; i < colored_unitigs.size(); i++) {
            if (i < (int64_t)colored_unitigs.size() - 1)  // Not last
                colored_unitigs[i].links.push_back(
                    colored_unitigs[i + 1].id);  // Chain unitigs
            else                                 // Last
                colored_unitigs[i].links =
                    U.links;  // The outgoing links of the original unitig
        }

        return colored_unitigs;
    }

    string get_unitig_string(vector<DBG::Node>& nodes, const DBG& dbg) {
        if (nodes.size() == 0) return "";
        string first_kmer = dbg.get_node_label(nodes[0]);
        string rest;
        for (int64_t i = 1; i < nodes.size(); i++)
            rest += dbg.incoming_character(nodes[i]);
        return first_kmer + rest;
    }

    void write_colorset(int64_t unitig_id, vector<uint32_t>& colorset,
                        ostream& colorsets_out) {
        colorsets_out << unitig_id;
        for (int64_t color : colorset) colorsets_out << " " << color;
        colorsets_out << "\n";
    }

    void write_unitig(vector<DBG::Node>& nodes, int64_t unitig_id, const DBG& dbg,
                      ostream& fasta_out, ostream& gfa_out) {
        string unitig_string = get_unitig_string(nodes, dbg);
        fasta_out << ">" << unitig_id << "\n" << unitig_string << "\n";
        gfa_out << "S\t" << unitig_id << "\t" << unitig_string << "\n";
    }

    void write_linkage(int64_t from_unitig, vector<int64_t> to_unitigs, ostream& gfa_out,
                       int64_t k) {
        for (int64_t to_unitig : to_unitigs) {
            // Print overlap with k-1 characters
            gfa_out << "L\t" << from_unitig << "\t+\t" << to_unitig << "\t+\t"
                    << (k - 1) << "M\n";
        }
    }

   public:
    // Writes the unitigs to the outputstream, one unitig per line
    // If split_by_colorset_runs == true, also splits the unitigs to maximal
    // runs of nodes that have the same colorset, and writes the color sets to
    // colorsets_out. If split_by_colorset_runs == false, then colorsets_out is
    // unused. If min_colors > 0, extracts maximal unitigs where each node has
    // at least >=min_colors colors. The unitigs are written in fasta-format.
    // The colorsets are written one per line, in the same order as unitigs, in
    // a space-separated format "id c1 c2 ..", where id is the fasta header of
    // the colorset, and c1 c2... are the colors.
    void extract_unitigs(const DBG& dbg, const Coloring& coloring, ostream& unitigs_out,
                         bool split_by_colorset_runs, ostream& colorsets_out,
                         ostream& gfa_out, int64_t min_colors = 0) {

        gfa_out << "H"
                << "\t"
                << "VN:Z:1.0"
                << "\n";  // Header

        unordered_map<DBG::Node, bool> visited;

        int64_t total_kmers = dbg.number_of_kmers();
        Progress_printer pp(total_kmers, 100);

        for(DBG::Node v : dbg.all_nodes()){
            if (!visited[v.id]) {
                Unitig U = get_unitig_containing_node(v, dbg, visited);
                for (int64_t i = 0; i < U.nodes.size(); i++)
                    pp.job_done();  // Record progress
                if (min_colors > 0) {
                    for (Colored_Unitig& CU :
                         split_by_colorset_size(min_colors, U, coloring)) {
                        write_unitig(CU.nodes, CU.id, dbg, unitigs_out,
                                     gfa_out);
                        write_linkage(CU.id, CU.links, gfa_out, dbg.get_k());
                    }
                } else if (split_by_colorset_runs) {
                    for (Colored_Unitig& CU :
                         split_to_colorset_runs(U, coloring)) {
                        write_unitig(CU.nodes, CU.id, dbg, unitigs_out,
                                     gfa_out);
                        write_linkage(CU.id, CU.links, gfa_out, dbg.get_k());
                        write_colorset(CU.id, CU.colorset, colorsets_out);
                    }
                } else {
                    write_unitig(U.nodes, U.id, dbg, unitigs_out, gfa_out);
                    write_linkage(U.id, U.links, gfa_out, dbg.get_k());
                }
            }
        }
    }
};
