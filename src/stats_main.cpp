#include "Themisto.hh"
#include "input_reading.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"
#include "extract_unitigs.hh"

using namespace std;

int stats_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Extract unitigs out of the Themisto index.");

    options.add_options()
        ("i,index-prefix", "The index prefix that was given to the build command.", cxxopts::value<string>())
        ("unitigs", "Also compute statistics on unitigs. This takes a while and requires the temporary directory to be set.", cxxopts::value<bool>()->default_value("false"))
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    string index_dbg_file = opts["index-prefix"].as<string>() + ".tdbg";
    string index_color_file = opts["index-prefix"].as<string>() + ".tcolors";
    bool do_unitigs = opts["unitigs"].as<bool>();
    check_readable(index_dbg_file);
    check_readable(index_color_file);

    if(do_unitigs) get_temp_file_manager().set_dir(opts["temp-dir"].as<string>());

    Themisto themisto;

    write_log("Loading the index", LogLevel::MAJOR);    
    themisto.load_boss(index_dbg_file);
    themisto.load_colors(index_color_file);

    write_log("Computing index statistics", LogLevel::MAJOR);
    vector<bool> dummy_marks = themisto.boss.get_dummy_node_marks();

    LL total_nodes = themisto.boss.number_of_nodes();
    LL dummy_nodes = 0;
    for(bool b : dummy_marks) if(b) dummy_nodes++;

    LL total_edges = themisto.boss.number_of_edges();
    LL dummy_edges = 0;
    for(LL v = 0; v < themisto.boss.number_of_nodes(); v++){
        if(dummy_marks[v]) dummy_edges += themisto.boss.outdegree(v);
    }

    cout << "Node length k: " << themisto.boss.get_k() << endl;
    cout << "Edge length k+1: " << themisto.boss.get_k() + 1 << endl;
    cout << "Node count: " << total_nodes - dummy_nodes << endl;
    cout << "Node count (including technical dummy nodes): " << total_nodes << endl;
    cout << "Edge count: " << total_edges - dummy_edges << endl;
    cout << "Edge count (including technical dummy edges): " << total_edges << endl;
    cout << "Number of colors: " << themisto.coloring.n_colors << endl;
    cout << "Number of distinct color sets: " << themisto.coloring.get_number_of_distinct_colorsets() << endl;
    cout << "Sum of sizes of all distinct color sets: " << themisto.coloring.get_color_set_concatenation_length() << endl;

    if(do_unitigs){
        write_log("Extracting unitigs (this could take a while)", LogLevel::MAJOR);
        UnitigExtractor UE;
        string unitigs_file = get_temp_file_manager().create_filename("unitigs-");
        throwing_ofstream unitigs_out(unitigs_file);
        NullStream null_stream;
        UE.extract_unitigs(themisto, unitigs_out.stream, false, null_stream);
        unitigs_out.close();
        
        LL unitig_count = 0;
        Sequence_Reader_Buffered sr(unitigs_file, FASTA_MODE);
        LL min_unitig_len = 1e18;
        LL max_unitig_len = 0;
        LL unitig_len_sum = 0;
        while(true){
            LL len = sr.get_next_read_to_buffer();
            if(len == 0) break;
            unitig_count++;
            min_unitig_len = min(min_unitig_len, len);
            max_unitig_len = max(max_unitig_len, len);
            unitig_len_sum += len;
        }

        cout << "Min unitig length: " << min_unitig_len << endl;    
        cout << "Max unitig length: " << max_unitig_len << endl;
        cout << "Avg unitig length: " << (double)unitig_len_sum / unitig_count << endl;

    }

    return 0;
    
}
