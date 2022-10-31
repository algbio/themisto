#include "coloring/Coloring.hh"
#include "globals.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"

using namespace sbwt;
using namespace std;

int color_set_diagnostics_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Prints stuff avoid the coloring data structure. For developers.");

    options.add_options()
        ("i,index-prefix", "The index prefix that was given to the build command.", cxxopts::value<string>())
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

    write_log("Loading the index", LogLevel::MAJOR);
    plain_matrix_sbwt_t SBWT;
    Coloring<> coloring;
    SBWT.load(index_dbg_file);
    coloring.load(index_color_file, SBWT);

    write_log("Computing stats", LogLevel::MAJOR);

    vector<int64_t> sum_of_bit_sizes(coloring.largest_color()+1); // 0..largest_color
    vector<int64_t> count_of_occurrences(coloring.largest_color()+1); // 0..largest_color

    for(const typename Coloring<>::colorset_type::view_t cs : coloring.get_all_distinct_color_sets()){
        sum_of_bit_sizes[cs.size()] += cs.size_in_bits();
        count_of_occurrences[cs.size()]++;
    }

    for(int64_t i = 0; i < sum_of_bit_sizes.size(); i++){
        cout << i << " " << count_of_occurrences[i] << " " << sum_of_bit_sizes[i] << endl;
    }

    return 0;
    
}
