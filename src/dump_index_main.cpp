
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "extract_unitigs.hh"
#include "DBG.hh"
#include "coloring/Coloring.hh"
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"
#include "new_extract_unitigs.hh"
#include "sbwt/cxxopts.hpp"

using namespace std;

int dump_index_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Dump index as maximal subunitigs where all k-mers have the same color set.");

    options.add_options()
        ("i,index-prefix", "The index prefix that was given to the build command.", cxxopts::value<string>())
        ("unitigs-out", "Output filename for the unitigs in FASTA format (optional).", cxxopts::value<string>())
        ("colors-out", "Output filename for the color sets (optional).", cxxopts::value<string>())
        ("metadata-out", "Output filename for the metadata (optional)", cxxopts::value<string>())
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

    optional<string> unitigs_outfile;
    try {
        unitigs_outfile = opts["unitigs-out"].as<string>();
    } catch(cxxopts::option_has_no_value_exception& e){
        // No file given -> ok.
    }

    optional<string> colors_outfile;
    try {
        colors_outfile = opts["colors-out"].as<string>();
    } catch(cxxopts::option_has_no_value_exception& e){
        // No file given -> ok.
    }

    optional<string> metadata_outfile;
    try {
        metadata_outfile = opts["metadata-out"].as<string>();
    } catch(cxxopts::option_has_no_value_exception& e){
        // No file given -> ok.
    }

    if(!unitigs_outfile.has_value() && !colors_outfile.has_value() && !metadata_outfile.has_value()){
        cerr << "ERROR: no output files given" << endl;
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
