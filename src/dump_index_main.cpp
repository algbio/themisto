
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "DBG.hh"
#include "coloring/Coloring.hh"
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"
#include "sbwt/cxxopts.hpp"
#include "dump_index.hh"

using namespace std;

int dump_index_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Dump index as maximal subunitigs where all k-mers have the same color set.");

    string unitig_file_suffix = ".unitigs.fa";
    string metadata_file_suffix = ".metadata.txt";
    string color_file_suffix = ".color_sets.txt";

    options.add_options()
        ("i,index-prefix", "The index prefix that was given to the build command.", cxxopts::value<string>())
        ("o,output-prefix", "Prefix for output filenames. Will write [prefix]" + unitig_file_suffix + ", [prefix]" + color_file_suffix + " and [prefix]" + metadata_file_suffix, cxxopts::value<string>())
        ("no-unitigs", "Do not dump unitigs.", cxxopts::value<bool>()->default_value("false"))
        ("no-color-sets", "Do not dump color sets.", cxxopts::value<bool>()->default_value("false"))
        ("no-metadata", "Do not dump metadata.", cxxopts::value<bool>()->default_value("false"))
        ("t, n-threads", "Number of parallel threads", cxxopts::value<int64_t>()->default_value("4"))
        ("v,verbose", "More verbose progress reporting into stderr.", cxxopts::value<bool>()->default_value("false"))
        ("silent", "Print as little as possible to stderr (only errors).", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    if(opts["verbose"].as<bool>() && opts["silent"].as<bool>())
        throw runtime_error("Can not give both --verbose and --silent");
    if(opts["verbose"].as<bool>()) set_log_level(LogLevel::MINOR);
    if(opts["silent"].as<bool>()) set_log_level(LogLevel::OFF);

    string out_prefix = opts["output-prefix"].as<string>();

    optional<string> unitigs_outfile;
    optional<string> colors_outfile;
    optional<string> metadata_outfile;

    if(!opts["no-unitigs"].as<bool>()) unitigs_outfile = out_prefix + unitig_file_suffix;
    if(!opts["no-color-sets"].as<bool>()) colors_outfile = out_prefix + color_file_suffix;
    if(!opts["no-metadata"].as<bool>()) metadata_outfile = out_prefix + metadata_file_suffix;

    if(!unitigs_outfile.has_value() && !colors_outfile.has_value() && !metadata_outfile.has_value()){
        cerr << "ERROR: all output was disabled" << endl;
        return 1;
    }

    int64_t n_threads = opts["n-threads"].as<int64_t>();
    string index_dbg_file = opts["index-prefix"].as<string>() + ".tdbg";
    string index_color_file = opts["index-prefix"].as<string>() + ".tcolors";

    check_readable(index_dbg_file);

    // Start

    write_log("Starting", LogLevel::MAJOR);
    write_log("Loading the index", LogLevel::MAJOR);

    sbwt::plain_matrix_sbwt_t SBWT;
    SBWT.load(index_dbg_file);

    std::variant<Coloring<SDSL_Variant_Color_Set>, Coloring<Roaring_Color_Set>> coloring;
    load_coloring(index_color_file, SBWT, coloring);

    write_log("Preparing the de Bruijn graph", LogLevel::MAJOR);
    DBG dbg(&SBWT);

    write_log("Dumping the index", LogLevel::MAJOR);
    std::visit([&](const auto& coloring){
        dump_index(n_threads, dbg, coloring, unitigs_outfile, colors_outfile, metadata_outfile);
    }, coloring);

    return 0;

}
