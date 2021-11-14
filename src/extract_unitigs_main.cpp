#include "Themisto.hh"
#include "input_reading.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"
#include "extract_unitigs.hh"

using namespace std;

int extract_unitigs_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Extract unitigs out of the Themisto index.");

    options.add_options()
        ("i,index-dir", "Location of the Themisto index.", cxxopts::value<string>())
        ("unitigs-out", "Output filename for the unitigs (outputted in FASTA format).", cxxopts::value<string>()->default_value(""))
        ("colors-out", "Output filename for the unitig colors. If this option is not given, the colors are not computed. Note that giving this option affects the unitigs written to unitigs-out: if a unitig has nodes with different color sets, the unitig is split into maximal segments of nodes that have equal color sets. The file format of the color file is as follows: there is one line for each unitig. The lines contain space-separated strings. The first string on a line is the FASTA header of a unitig, and the following strings on the line are the integer color labels of the colors of that unitig. The unitigs appear in the same order as in the FASTA file.", cxxopts::value<string>()->default_value(""))
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    string index_dir = opts["index-dir"].as<string>();
    check_true(index_dir != "", "Index directory not set");
    check_dir_exists(index_dir);

    string unitigs_outfile = opts["unitigs-out"].as<string>();
    string colors_outfile = opts["colors-out"].as<string>();
    bool do_colors = (colors_outfile != "");
    if(unitigs_outfile == "" && colors_outfile == ""){
        write_log("Error: no output files given");
        return 1;
    }
    
    // Prepare output streams

    NullStream null_stream;
    throwing_ofstream unitigs_ofstream;
    throwing_ofstream colors_ofstream;

    ostream* unitigs_out = &null_stream;
    ostream* colors_out = &null_stream;

    if(unitigs_outfile != ""){
        unitigs_ofstream.open(unitigs_outfile);
        unitigs_out = &(unitigs_ofstream.stream);
    }

    if(colors_outfile != ""){
        colors_ofstream.open(colors_outfile);
        colors_out = &(colors_ofstream.stream);
    }

    // Start

    write_log("Starting");
    write_log("Loading the index");    
    Themisto themisto;
    themisto.load_boss(index_dir + "/boss-");
    themisto.load_colors(index_dir + "/coloring-");

    write_log("Extracting unitigs");
    
    UnitigExtractor UE;
    UE.extract_unitigs(themisto, *unitigs_out, do_colors, *colors_out);

    return 0;
    
}
