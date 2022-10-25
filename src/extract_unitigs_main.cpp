#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"
#include "extract_unitigs.hh"
#include "DBG.hh"
#include "Coloring.hh"
#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"

using namespace std;

int extract_unitigs_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Extract unitigs out of the Themisto index.");

    options.add_options()
        ("i,index-prefix", "The index prefix that was given to the build command.", cxxopts::value<string>())
        ("fasta-out", "Output filename for the unitigs in FASTA format (optional).", cxxopts::value<string>()->default_value(""))
        ("gfa-out", "Output the unitig graph in GFA1 format (optional).", cxxopts::value<string>()->default_value(""))
        ("colors-out", "Output filename for the unitig colors (optional). If this option is not given, the colors are not computed. Note that giving this option affects the unitigs written to unitigs-out: if a unitig has nodes with different color sets, the unitig is split into maximal segments of nodes that have equal color sets. The file format of the color file is as follows: there is one line for each unitig. The lines contain space-separated strings. The first string on a line is the FASTA header of a unitig (without the '>'), and the following strings on the line are the integer color labels of the colors of that unitig. The unitigs appear in the same order as in the FASTA file.", cxxopts::value<string>()->default_value(""))
        ("min-colors", "Extract maximal unitigs with at least (>=) min-colors in each node. Can't be used with --colors-out. (optional)", cxxopts::value<LL>()->default_value("0"))
        ("v,verbose", "More verbose progress reporting into stderr.", cxxopts::value<bool>()->default_value("false"))
        ("silent", "Print as little as possible to stderr (only errors).", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    if(opts["verbose"].as<bool>() && opts["silent"].as<bool>())
        throw runtime_error("Can not give both --verbose and --silent");
    if(opts["verbose"].as<bool>()) set_log_level(LogLevel::MINOR);
    if(opts["silent"].as<bool>()) set_log_level(LogLevel::OFF);

    string unitigs_outfile = opts["fasta-out"].as<string>();
    string gfa_outfile = opts["gfa-out"].as<string>();
    string colors_outfile = opts["colors-out"].as<string>();
    LL min_colors = opts["min-colors"].as<LL>();
    bool do_colors = (colors_outfile != "") || min_colors > 0;
    string index_dbg_file = opts["index-prefix"].as<string>() + ".tdbg";
    string index_color_file = opts["index-prefix"].as<string>() + ".tcolors";

    check_readable(index_dbg_file);
    if(do_colors) check_readable(index_color_file);

    if(unitigs_outfile == "" && colors_outfile == "" && gfa_outfile == ""){
        throw std::runtime_error("Error: no output files given");
    }
    if (colors_outfile != "" && min_colors != 0) { // Colors may not make sense for --min-colors
    throw std::runtime_error("Error: Colorset extraction not allowed when --min-colors is given");
    }

    // Prepare output streams

    sbwt::SeqIO::NullStream null_stream;
    throwing_ofstream unitigs_ofstream;
    throwing_ofstream gfa_ofstream;
    throwing_ofstream colors_ofstream;

    ostream* unitigs_out = &null_stream;
    ostream* gfa_out = &null_stream;
    ostream* colors_out = &null_stream;

    if(unitigs_outfile != ""){
        unitigs_ofstream.open(unitigs_outfile);
        unitigs_out = &(unitigs_ofstream.stream);
    }

    if(gfa_outfile != ""){
        gfa_ofstream.open(gfa_outfile);
        gfa_out = &(gfa_ofstream.stream);
    }

    if(colors_outfile != ""){
        colors_ofstream.open(colors_outfile);
        colors_out = &(colors_ofstream.stream);
    }

    // Start

    write_log("Starting", LogLevel::MAJOR);
    write_log("Loading the index", LogLevel::MAJOR);

    sbwt::plain_matrix_sbwt_t SBWT;
    SBWT.load(index_dbg_file);

    std::variant<Coloring<Color_Set>, Coloring<Roaring_Color_Set>, Coloring<Bit_Magic_Color_Set>> coloring;
    if(do_colors){
        // Load whichever coloring data structure type is stored on disk
        load_coloring(index_color_file, SBWT, coloring);
    }

    DBG dbg(&SBWT);

    write_log("Extracting unitigs", LogLevel::MAJOR);

    if(std::holds_alternative<Coloring<Color_Set>>(coloring)){
        UnitigExtractor<Coloring<Color_Set>> UE;
        UE.extract_unitigs(dbg, std::get<Coloring<Color_Set>>(coloring), *unitigs_out, do_colors, *colors_out, *gfa_out, min_colors);
    }
    if(std::holds_alternative<Coloring<Roaring_Color_Set>>(coloring)){
        UnitigExtractor<Coloring<Roaring_Color_Set>> UE;
        UE.extract_unitigs(dbg, std::get<Coloring<Roaring_Color_Set>>(coloring), *unitigs_out, do_colors, *colors_out, *gfa_out, min_colors);
    }
    if(std::holds_alternative<Coloring<Bit_Magic_Color_Set>>(coloring)){
        UnitigExtractor<Coloring<Bit_Magic_Color_Set>> UE;
        UE.extract_unitigs(dbg, std::get<Coloring<Bit_Magic_Color_Set>>(coloring), *unitigs_out, do_colors, *colors_out, *gfa_out, min_colors);
    }


    return 0;

}
