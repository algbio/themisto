#include "sbwt/SBWT.hh"
#include "sbwt/globals.hh"
#include "globals.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <variant>
#include "version.h"
#include "cxxopts.hpp"
#include "DBG.hh"
#include "coloring/Coloring.hh"

using namespace sbwt;
using namespace std;

int64_t count_DBG_edges(const DBG& dbg){
    int64_t ans = 0;
    for(DBG::Node v : dbg.all_nodes()){
        ans += dbg.outdegree(v);
    }
    return ans;
}

string human_readable_bytes(int64_t bytes){
    stringstream ss;
    ss << std::fixed << std::setprecision(3); // 3 digits after decimal point
    if(bytes < (1<<10)) ss << bytes << " bytes";
    else if(bytes < (1<<20)) ss << (double)bytes / (1 << 10) << " kB";
    else if(bytes < (1<<30)) ss << (double)bytes / (1 << 20) << " MB";
    else ss << (double)bytes / (1 << 30) << " GB";
    return ss.str();
}

int stats_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Compute statistics from a Themisto index.");

    options.add_options()
        ("i,index-prefix", "The index prefix that was given to the build command.", cxxopts::value<string>())
        ("space-breakdown", "Also give a space breakdown for the components of the index.", cxxopts::value<bool>()->default_value("false"))
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
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
    bool do_unitigs = opts["unitigs"].as<bool>();
    bool space_breakdown = opts["space-breakdown"].as<bool>();

    if(opts["verbose"].as<bool>() && opts["silent"].as<bool>())
        throw runtime_error("Can not give both --verbose and --silent");
    if(opts["verbose"].as<bool>()) set_log_level(LogLevel::MINOR);
    if(opts["silent"].as<bool>()) set_log_level(LogLevel::OFF);

    check_readable(index_dbg_file);
    check_readable(index_color_file);

    if(do_unitigs) get_temp_file_manager().set_dir(opts["temp-dir"].as<string>());

    write_log("Loading the SBWT", LogLevel::MAJOR);

    plain_matrix_sbwt_t SBWT;
    SBWT.load(index_dbg_file);

    cout << "Node length k: " << SBWT.get_k() << endl;
    cout << "Number of k-mers: " << SBWT.number_of_kmers() << endl;
    cout << "Number of subsets in the SBWT data structure: " << SBWT.number_of_subsets() << endl;

    std::variant<Coloring<SDSL_Variant_Color_Set>, Coloring<Roaring_Color_Set>> coloring;
    load_coloring(index_color_file, SBWT, coloring);

    if(std::holds_alternative<Coloring<SDSL_Variant_Color_Set>>(coloring))
        write_log("sdsl coloring structure loaded", LogLevel::MAJOR);
    if(std::holds_alternative<Coloring<Roaring_Color_Set>>(coloring))
        write_log("roaring coloring structure loaded", LogLevel::MAJOR);

    // Helper functions to be able to call member functions of coloring with std::visit.
    // This cleans up the code so that we don't have the branch where we check which
    // variant we currently have.
    auto call_largest_color = [](auto& obj) { return obj.largest_color(); };
    auto call_number_of_distinct_color_sets = [](auto& obj) { return obj.number_of_distinct_color_sets(); };
    auto call_sum_of_all_distinct_color_set_lengths = [](auto& obj) { return obj.sum_of_all_distinct_color_set_lengths(); };
    auto call_space_breakdown = [](auto& obj) { return obj.space_breakdown(); };

    cout << "Color id range: 0.." << std::visit(call_largest_color, coloring) << endl;
    cout << "Number of distinct color sets: " << std::visit(call_number_of_distinct_color_sets, coloring) << endl;
    cout << "Sum of sizes of all distinct color sets: " << std::visit(call_sum_of_all_distinct_color_set_lengths, coloring) << endl;

    write_log("Computing more statistics...", LogLevel::MAJOR);
    DBG dbg(&SBWT);
    cout << "De Bruijn graph edge count: " << count_DBG_edges(dbg) << endl;

    if(space_breakdown){
        cout << "== Space breakdown of the coloring structure ==" << endl;
        for(auto [component, space] : std::visit(call_space_breakdown, coloring)){
            cout << component << ": " << human_readable_bytes(space) << endl;
        }
        cout << "== Space taken for the de Bruijn graph ==" << endl;
        seq_io::NullStream ns;
        int64_t bytes = SBWT.serialize(ns);
        cout << "SBWT: " << human_readable_bytes(bytes) << endl;
        cout << "==" << endl;

    }

    return 0;

}
