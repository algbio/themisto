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

template<typename coloring_t> 
void dump_distinct_color_sets_to_binary(const coloring_t& coloring, const string& outfile){

    int64_t n = coloring.number_of_distinct_color_sets();

    Buffered_ofstream<> out(outfile);

    vector<int64_t> color_buf;
    sbwt::Progress_printer pp(n, 100);
    for(int64_t color_set_id = 0; color_set_id < n; color_set_id++){
        color_buf.clear();
        typename coloring_t::colorset_view_type cs = coloring.get_color_set_by_color_set_id(color_set_id);
        cs.push_colors_to_vector(color_buf);
        if(cs.size() > UINT32_MAX){
            throw std::runtime_error("Error: color set has more than 2^32 elements");
        }
        uint32_t size = color_buf.size();
        out.write((char*)&size, sizeof(uint32_t));

        for(int64_t i = 0; i < color_buf.size(); i++){
            if(color_buf[i] > UINT32_MAX){
                throw std::runtime_error("Error: color id is larger than 2^32");
            }
            uint32_t color = color_buf[i];
            out.write((char*)&color, sizeof(uint32_t));
        }
    }

}

int dump_distinct_color_sets_to_binary_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Dumps distinct color sets in binary format where each color set is encoded with first a 32-bit integer describing the number of colors in the set, and then the colors as 32-bit integers.");

    options.add_options()
        ("i,index-prefix", "The index prefix that was given to the build command.", cxxopts::value<string>())
        ("o,output-file", "The output file for the dump.", cxxopts::value<string>())
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

    string index_dbg_file = opts["index-prefix"].as<string>() + ".tdbg";
    string index_color_file = opts["index-prefix"].as<string>() + ".tcolors";
    string outfile = opts["output-file"].as<string>();
    
    if(opts["verbose"].as<bool>()) set_log_level(LogLevel::MINOR);
    if(opts["silent"].as<bool>()) set_log_level(LogLevel::OFF);

    check_readable(index_dbg_file);
    check_readable(index_color_file);

    write_log("Loading the index", LogLevel::MAJOR);

    plain_matrix_sbwt_t SBWT;
    SBWT.load(index_dbg_file);
    
    std::variant<Coloring<SDSL_Variant_Color_Set>, Coloring<Roaring_Color_Set>> coloring;
    load_coloring(index_color_file, SBWT, coloring);

    auto call_dump_colors = [&](const auto& obj){
        dump_distinct_color_sets_to_binary(obj, outfile);
    };

    write_log("Dumping distinct color sets to " + outfile, LogLevel::MAJOR);
    std::visit(call_dump_colors, coloring);

    return 0;

}
