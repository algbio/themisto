#include "sbwt/SBWT.hh"
#include "sbwt/globals.hh"
#include "globals.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include <sstream>
#include <iomanip>
#include "version.h"
#include "cxxopts.hpp"
#include "extract_unitigs.hh"
#include "coloring/Coloring.hh"
#include "DBG.hh"

using namespace sbwt;
using namespace std;

// Assumes the row has been initialized with enough space for all colors
template<typename coloring_t> 
void print_color_set_as_bitmap(int64_t node_id, const coloring_t& coloring, Buffered_ofstream<>& out){

    char space = ' '; out.write(&space, 1); // Write a space

    vector<char> row(coloring.largest_color() + 1, '0'); // ASCII characters '0' and '1'

    // Set the colors
    for(int64_t color : coloring.get_color_set_of_node(node_id).get_colors_as_vector()) 
        row[color] = '1';

    // Write out
    out.write(row.data(), row.size());
}

// Assumes the row has been initialized with enough space for all colors
template<typename coloring_t> 
void print_color_set_as_integers(int64_t node_id, const coloring_t& coloring, Buffered_ofstream<>& out){
    char string_buf[32]; // Enough space to represent a 64-bit integer in ascii
    for(int64_t color : coloring.get_color_set_of_node(node_id).get_colors_as_vector()){
        int64_t len = fast_int_to_string(color, string_buf);

        char space = ' '; out.write(&space, 1); 
        out.write(string_buf, len);
    }
}

template<typename coloring_t> 
void dump_colors(const DBG& dbg, const coloring_t& coloring, string outfile, bool sparse){

    int64_t n_colors = coloring.largest_color() + 1;

    Buffered_ofstream<> out(outfile);
    vector<bool> matrix_row(n_colors, 0);
    vector<char> matrix_row_as_ASCII(n_colors, '0');

    sbwt::Progress_printer pp(dbg.number_of_kmers(), 100);
    for(DBG::Node v : dbg.all_nodes()){
        // Extract the k-mer
        string kmer = dbg.get_node_label(v);

        // Write the k-mer and a space
        out.write(kmer.c_str(), kmer.size());

        if(sparse) print_color_set_as_integers(v.id, coloring, out);
        else print_color_set_as_bitmap(v.id, coloring, out);

        char newline = '\n'; out.write(&newline, 1); // Write a line break

        pp.job_done();
    }

}

int dump_color_matrix_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "This command prints a file where each line corresponds to a k-mer in the index. The line starts with the k-mer, followed by space, followed by the color set of that k-mer. If `--sparse` is given, the color set is printed as a space-separated list of integers. Otherwise, the color set is printed as a string of zeroes and ones such that the i-th character is '1' iff color i is present in the color set.");

    options.add_options()
        ("i,index-prefix", "The index prefix that was given to the build command.", cxxopts::value<string>())
        ("o,output-file", "The output file for the dump.", cxxopts::value<string>())
        ("v,verbose", "More verbose progress reporting into stderr.", cxxopts::value<bool>()->default_value("false"))
        ("silent", "Print as little as possible to stderr (only errors).", cxxopts::value<bool>()->default_value("false"))
        ("sparse", "Print only the indices of non-zero entries.", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    string index_dbg_file = opts["index-prefix"].as<string>() + ".tdbg";
    string index_color_file = opts["index-prefix"].as<string>() + ".tcolors";
    string outfile = opts["output-file"].as<string>();
    bool sparse = opts["sparse"].as<bool>();

    if(opts["verbose"].as<bool>()) set_log_level(LogLevel::MINOR);
    if(opts["silent"].as<bool>()) set_log_level(LogLevel::OFF);

    check_readable(index_dbg_file);
    check_readable(index_color_file);

    write_log("Loading the index", LogLevel::MAJOR);

    plain_matrix_sbwt_t SBWT;
    SBWT.load(index_dbg_file);
    DBG dbg(&SBWT);

    std::variant<Coloring<SDSL_Variant_Color_Set>, Coloring<Roaring_Color_Set>> coloring;
    load_coloring(index_color_file, SBWT, coloring);

    auto call_dump_colors = [&](const auto& obj){
        dump_colors(dbg, obj, outfile, sparse);
    };

    write_log("Dumping color matrix to " + outfile, LogLevel::MAJOR);
    std::visit(call_dump_colors, coloring);

    return 0;

}
