#pragma once

#include <iostream>
#include "globals.hh"
#include "Themisto.hh"

using namespace std;

class UnitigExtractor{

private:

    struct Unitig{
        vector<LL> nodes;
        vector<LL> links; // Outgoing edges from this unitig, as a list of unitig ids
        LL id;
    };

    struct Colored_Unitig{
        vector<LL> nodes;
        vector<LL> colorset;
        vector<LL> links; // Outgoing edges from this unitig, as a list of unitig ids
        LL id;
    };

    // Assumes that dummy nodes are always visited
    Unitig get_node_unitig_containing(LL v, const BOSS<sdsl::bit_vector>& boss, vector<bool>& visited){
        
        vector<LL> backward, forward;

        visited[v] = true;

        // Walk backward until (indegree >= 2), or (predecessor is visited or has outdegree >= 2)
        LL u = v;
        while(true){
            if(boss.indegree(u) >= 2) break;
            LL pred = boss.edge_source(boss.inedge_range(u).first); // predecessor
            if(pred == -1) break;
            if(visited[pred]) break;
            if(boss.outdegree(pred) >= 2) break;
            u = pred;
            visited[u] = true;
            backward.push_back(u);
        }

        // Walk forward until (outdegree >= 2) or (successor is visited or has indegree >= 2)
        u = v;
        while(true){
            if(boss.outdegree(u) >= 2) break;
            LL succ = boss.walk(u, boss.outlabels_at(boss.outlabel_range(u).first));
            if(succ == -1) break;
            if(visited[succ]) break;
            if(boss.indegree(succ) >= 2) break;
            u = succ;
            visited[u] = true;
            forward.push_back(u);
        }

        // Collect the unitig path
        Unitig U;
        for(LL i = (LL)backward.size()-1; i >= 0; i--) U.nodes.push_back(backward[i]);
        U.nodes.push_back(v);
        for(LL w : forward) U.nodes.push_back(w);
        
        // Compute links
        pair<LL,LL> outlabel_range = boss.outlabel_range(U.nodes.back());
        for(LL i = outlabel_range.first; i <= outlabel_range.second; i++){
            LL successor = boss.walk(U.nodes.back(), boss.outlabels_at(i));
            assert(successor != -1); // Must exist
            U.links.push_back(successor);
        }

        U.id = U.nodes[0];
        return U;
    }

    vector<Colored_Unitig> split_to_colorset_runs(Unitig& U, Themisto& themisto){
        vector<Colored_Unitig> colored_unitigs;
        auto get_id = [&](LL node){return themisto.coloring.get_colorset_id(node, themisto.boss);};
        vector<LL> colorset_buffer(themisto.coloring.n_colors);
        for(LL run_start = 0; run_start < U.nodes.size(); run_start++){
            LL run_end = run_start;
            while(run_end < U.nodes.size()-1 && get_id(U.nodes[run_start]) == get_id(U.nodes[run_end+1])){
                run_end++;
            }

            Colored_Unitig colored;
            for(LL i = run_start; i <= run_end; i++) 
                colored.nodes.push_back(U.nodes[i]);
            colored.colorset = themisto.coloring.get_colorvec(U.nodes[run_start], themisto.boss);
            colored.id = U.nodes[run_start];

            colored_unitigs.push_back(colored); // Links are added later
            
            run_start = run_end;
        }

        // Linkage
        for(LL i = 0; i < colored_unitigs.size(); i++){
            if(i < (LL)colored_unitigs.size()-1) // Not last
                colored_unitigs[i].links.push_back(colored_unitigs[i].id); // Chain unitigs
            else // Last
                colored_unitigs[i].links = U.links; // The outgoing links of the original unitig
        }

        return colored_unitigs;
    }

    string get_unitig_string(vector<LL>& nodes, Themisto& themisto){
        if(nodes.size() == 0) return "";
        string first_kmer = themisto.boss.get_node_label(nodes[0]);
        string rest;
        for(LL i = 1; i < nodes.size(); i++)
            rest += themisto.boss.incoming_character(nodes[i]);
        return first_kmer + rest;
    }

    void write_colorset(LL unitig_id, vector<LL>& colorset, ostream& colorsets_out){
        colorsets_out << unitig_id;
        for(LL color : colorset) colorsets_out << " " << color;
        colorsets_out << "\n";
    }

    void write_unitig(vector<LL>& nodes, LL unitig_id, Themisto& themisto, ostream& fasta_out, ostream& gfa_out){
        string unitig_string = get_unitig_string(nodes, themisto);
        fasta_out << ">" << unitig_id<< "\n" << unitig_string << "\n";
        gfa_out << "S\t" << unitig_id << "\t" << unitig_string << "\n";
    }

    void write_linkage(LL from_unitig, vector<LL> to_unitigs, ostream& gfa_out, LL k){
        for(LL to_unitig : to_unitigs){
            // Print overlap with k-1 characters
            gfa_out << "L\t" << from_unitig << "\t+\t" << to_unitig << "\t+\t " << (k-1) << "M\n";
        }
    }

public:

    // Writes the unitigs to the outputstream, one unitig per line
    // If split_by_colorset_runs == true, also splits the unitigs to maximal runs of nodes that have
    // the same colorset, and writes the color sets to colorsets_out. If split_by_colorset_runs == false,
    // then colorsets_out is unused.
    // The unitigs are written in fasta-format.
    // The colorsets are written one per line, in the same order as unitigs, in a space-separated format "id c1 c2 ..",
    // where id is the fasta header of the colorset, and c1 c2... are the colors.
    void extract_unitigs(Themisto& themisto, ostream& unitigs_out, bool split_by_colorset_runs, ostream& colorsets_out, ostream& gfa_out){
        BOSS<sdsl::bit_vector>& boss = themisto.boss;

        gfa_out << "H" << "\t" << "VN:Z:1.0" << "\n"; // Header

        vector<bool> visited = boss.get_dummy_node_marks();

        LL non_dummies = 0;
        for(LL i = 0; i < visited.size(); i++) if(!visited[i]) non_dummies++;
        Progress_printer pp(non_dummies, 100);

        for(LL v = 0; v < themisto.boss.number_of_nodes(); v++){
            if(!visited[v]){
                Unitig U = get_node_unitig_containing(v, boss, visited);
                for(LL i = 0; i < U.nodes.size(); i++) pp.job_done(); // Record progress
                if(split_by_colorset_runs){
                    for(Colored_Unitig& CU : split_to_colorset_runs(U, themisto)){
                        write_unitig(CU.nodes, CU.id, themisto, unitigs_out, gfa_out);
                        write_linkage(CU.id, CU.links, gfa_out, boss.get_k());
                        write_colorset(CU.id, CU.colorset, colorsets_out);
                    }
                } else{
                    write_unitig(U.nodes, U.id, themisto, unitigs_out, gfa_out);
                    write_linkage(U.id, U.links, gfa_out, boss.get_k());
                }
            }
        }
    }

};
