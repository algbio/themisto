#include "coloring/Coloring.hh"
#include "globals.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"

using namespace sbwt;
using namespace std;

int make_d_equal_1_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Makes d = 1. For developers.");

    options.add_options()
        ("i", "The index prefix that was given to the build command.", cxxopts::value<string>())
        ("o", "The index prefix for the output.", cxxopts::value<string>())
        ("t", "Number of threads.", cxxopts::value<int64_t>()->default_value("1"))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    string input_dbg_file = opts["i"].as<string>() + ".tdbg";
    string input_color_file = opts["i"].as<string>() + ".tcolors";

    string output_dbg_file = opts["o"].as<string>() + ".tdbg";
    string output_color_file = opts["o"].as<string>() + ".tcolors";

    int64_t n_threads = opts["t"].as<int64_t>();

    write_log("Loading the index", LogLevel::MAJOR);
    plain_matrix_sbwt_t SBWT;
    Coloring<> coloring;

    cerr << "Loading SBWT" << endl;
    SBWT.load(input_dbg_file);
    cerr << "Loading coloring" << endl;
    coloring.load(input_color_file, SBWT);

    cerr << "Building backward traversal support" << endl;
    SBWT_backward_traversal_support backward_support(&SBWT);

    cerr << "Making d = 1" << endl;
    coloring.add_all_node_id_to_color_set_id_pointers(SBWT, backward_support, n_threads);

    write_log("Saving the updated index", LogLevel::MAJOR);

    SBWT.serialize(output_dbg_file);
    coloring.serialize(output_color_file);

    write_log("Done", LogLevel::MAJOR);

    return 0;
    
}