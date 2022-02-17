#include "Themisto.hh"
#include "input_reading.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"

using namespace std;

vector<bool> mark_first_incoming_edges(BOSS<sdsl::bit_vector>& boss){
    // Mark the colex-smallest edge to each node.
    LL node_id = -1; // Wheeler rank
    LL edge_count_to_current_node = 0;
    LL count_of_current_inlabel = 0;
    char current_label = '\0';
    vector<bool> marks(boss.number_of_edges());
    Progress_printer pp(boss.indegs_size(), 100);
    for(LL i = 0; i < boss.indegs_size(); i++){
        if(boss.indegs_at(i) == 1){
            node_id++;
            edge_count_to_current_node = 0;
        }
        else{
            // Process edge
            char label = boss.incoming_character(node_id);
            if(label != current_label) count_of_current_inlabel = 0;
            count_of_current_inlabel++;
            
            edge_count_to_current_node++;

            if(edge_count_to_current_node == 1){
                // Find the edge in the outlabels string and mark.
                marks[boss.outlabels_select(count_of_current_inlabel, label)] = 1;
            }
        }
        pp.job_done();
    }
    return marks;
}

void dump_nodeboss(BOSS<sdsl::bit_vector>& boss, string out_prefix){

}

int main(int argc, char** argv){

    set_log_level(LogLevel::MINOR);

    // Legacy support: transform old option format --outfile --out-file
    string legacy_support_fix = "--out-file";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--outfile") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "This program prints a NodeBOSS representation of the index.");

    options.add_options()
        ("o,out-prefix", "Output file prefix.", cxxopts::value<string>())
        ("i,dbg-file", "The .tdbg file of an index.", cxxopts::value<string>())
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        cerr << "Usage examples: " << argv[0] << " -i in_prefix -o out_prefix --temp-dir temp" << endl;
        exit(1);
    }

    string out_prefix = opts["out-prefix"].as<string>();
    string index_dbg_file = opts["dbg-file"].as<string>();
    string temp_dir = opts["temp-dir"].as<string>();

    create_directory_if_does_not_exist(temp_dir);

    write_log("Starting", LogLevel::MAJOR);

    get_temp_file_manager().set_dir(temp_dir);

    write_log("Loading the DBG", LogLevel::MAJOR);
    Themisto themisto;
    themisto.load_boss(index_dbg_file);

    write_log("Marking edges to keep", LogLevel::MAJOR);
    vector<bool> marks = mark_first_incoming_edges(themisto.boss);

    return 0;
}
