#include "coloring/Coloring.hh"
#include "globals.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"
#include "pseudoalign.hh"

using namespace sbwt;
using namespace std;

int query_color_set_ids_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Prints the color set ids of the query strings to stdout. For developers.");

    options.add_options()
        ("i", "The index prefix that was given to the build command.", cxxopts::value<string>())
        ("q", "The query file (fasta or fastq, possibly gzipped)", cxxopts::value<string>())
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
    string query_file = opts["q"].as<string>();

    check_readable(query_file);

    write_log("Loading the index", LogLevel::MAJOR);
    plain_matrix_sbwt_t SBWT;
    Coloring<> coloring;

    cerr << "Loading SBWT" << endl;
    SBWT.load(input_dbg_file);
    cerr << "Loading coloring" << endl;
    coloring.load(input_color_file, SBWT);

    cerr << "Running queries" << endl;

    if(seq_io::figure_out_file_format(query_file).gzipped){
        seq_io::Reader<seq_io::Buffered_ifstream<seq_io::zstr::ifstream>> reader(query_file);
        print_query_color_set_ids(SBWT, coloring, reader);
    } else{
        seq_io::Reader<seq_io::Buffered_ifstream<std::ifstream>> reader(query_file);
        print_query_color_set_ids(SBWT, coloring, reader);
    }

    write_log("Done", LogLevel::MAJOR);

    return 0;
    
}