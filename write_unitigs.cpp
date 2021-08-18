#include "Themisto.hh"
#include "input_reading.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"

using namespace std;

char get_rc(char c){
    switch(c){
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return c;
    }
}

vector<string> read_lines(string filename){
    check_readable(filename);
    vector<string> lines;
    throwing_ifstream in(filename);
    string line;
    while(in.getline(line)){
        lines.push_back(line);
    }
    return lines;
}

struct Unitig{
    vector<LL> nodes;
    vector<char> node_labels;
};

// Assumes v is a node that starts a unitig, i.e. either
// - v has indegree 0, or
// - v has indegree 1 and the predecessor of v has indegree or outdegree >= 2
// Returns the longest path from v such that all nodes on the path (including v itself)
// have indegree and outdegree at most 1. May return an empty path.
Unitig get_node_unitig_from(LL v, const BOSS<sdsl::bit_vector>& boss){
    Unitig U;
    while(boss.indegree(v) <= 1 && boss.outdegree(v) <= 1){
        // v is a node in the unitig
        U.nodes.push_back(v);
        U.node_labels.push_back(boss.incoming_character(v));
        if(boss.outdegree(v) == 1)
            v = boss.walk(v, boss.node_outlabels(v)[0]);
        else break;
    }
    return U;
}

void split_to_colorset_runs_and_write(Unitig& U, Themisto& themisto, throwing_ofstream& out){
    if(U.nodes.size() == 0) return;
    auto get_id = [&](LL node){return themisto.coloring.get_colorset_id(node, themisto.boss);};
    vector<LL> colorset_buffer(themisto.coloring.n_colors);
    for(LL run_start = 0; run_start < U.nodes.size(); run_start++){
        LL run_end = run_start;
        while(run_end < U.nodes.size()-1 && get_id(U.nodes[run_start]) == get_id(U.nodes[run_end+1])){
            run_end++;
        }

        if(run_end - run_start + 1 >= themisto.boss.get_k()){
            // Write the unitig section and the colors
            LL n_colors = themisto.coloring.get_colorset_to_buffer(U.nodes[run_start], themisto.boss, colorset_buffer);
            if(n_colors > 0){
                for(LL i = run_start; i <= run_end; i++) out << U.node_labels[i];
                out << " ";
                for(LL i = 0; i < n_colors; i++) out << colorset_buffer[i] << " ";
                out << "\n";
            }
        }
        run_start = run_end;
    }

}

int main2(int argc, char** argv){

    // Legacy support: transform old option format --outfile --out-file
    string legacy_support_fix = "--out-file";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--outfile") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "");

    options.add_options()
        ("i,index-dir", "Directory where the index will be built. Always required, directory must exist before running.", cxxopts::value<string>())
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    // todo: positional argument

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    string index_dir = opts["index-dir"].as<string>();
    check_true(index_dir != "", "Index directory not set");
    check_dir_exists(index_dir);

    string outfile = opts["out-file"].as<string>();
    check_writable(outfile);

    write_log("Starting");

    //temp_file_manager.set_dir(""); // There should not be any temp files

    write_log("Loading the index");    
    Themisto themisto;
    themisto.load_boss(index_dir + "/boss-");
    themisto.load_colors(index_dir + "/coloring-");

    write_log("Writing the unitigs");
    throwing_ofstream out(outfile);

    // A unitig is a maximal non-branching path of edges, i.e.
    // a path of edges such that every node in the path except the first
    // and the last node have in-degree and out-degree exactly 1, and the
    // endpoints have in-edgee our outdegree 0 or >= 2

    // The unitigs partition the edges: each edge is part of exactly one unitig.

    // The nodes of the unitig are all the nodes in the path except the endpoints

    // We split the unitigs into runs with equal color set

    BOSS<sdsl::bit_vector>& boss = themisto.boss;
    Progress_printer pp(themisto.boss.number_of_nodes(), 100);
    for(LL v = 0; v < themisto.boss.number_of_nodes(); v++){
        Unitig U;
        if(boss.indegree(v) == 0){
            U = get_node_unitig_from(v, boss);
        } else if(boss.indegree(v) == 1){
            LL u = boss.edge_source(boss.inedge_range(v).first); // predecessor
            if(boss.indegree(u) >= 2 || boss.outdegree(u) >= 2){
                U = get_node_unitig_from(v, boss);
            }
        }
        split_to_colorset_runs_and_write(U, themisto, out);
        pp.job_done();
    }

    write_log("Finished");

    return 0;
}

int main(int argc, char** argv){
    write_log("Themisto-" + std::string(THEMISTO_BUILD_VERSION));
    write_log("Maximum k-mer length (size of the de Bruijn graph node labels): " + std::to_string(KMER_MAX_LENGTH-1));
    try{
        return main2(argc, argv);
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    }
}
