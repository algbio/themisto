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

    cxxopts::Options options(argv[0], "Extract unitigs out of the Themisto index.");

    options.add_options()
        ("i,index-prefix", "The index prefix that was given to the build command.", cxxopts::value<string>())
        ("unitigs", "Also compute statistics on unitigs. This takes a while and requires the temporary directory to be set.", cxxopts::value<bool>()->default_value("false"))
        ("space-breakdown", "Also give a space breakdown for the components of the index.", cxxopts::value<bool>()->default_value("false"))
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
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

    write_log("Loading the index", LogLevel::MAJOR);

    plain_matrix_sbwt_t SBWT;
    SBWT.load(index_dbg_file);

    std::variant<Coloring<Color_Set>, Coloring<Roaring_Color_Set>, Coloring<Bit_Magic_Color_Set>> coloring;
    load_coloring(index_color_file, SBWT, coloring);

    if(std::holds_alternative<Coloring<Color_Set>>(coloring))
        write_log("sdsl coloring structure loaded", LogLevel::MAJOR);
    if(std::holds_alternative<Coloring<Roaring_Color_Set>>(coloring))
        write_log("roaring coloring structure loaded", LogLevel::MAJOR);
    if(std::holds_alternative<Coloring<Bit_Magic_Color_Set>>(coloring))
        write_log("BitMagic coloring structure loaded", LogLevel::MAJOR);

    // Helper functions to be able to call member functions of coloring with std::visit.
    // This cleans up the code so that we don't have the branch where we check which
    // variant we currently have.
    auto call_largest_color = [](auto& obj) { return obj.largest_color(); };
    auto call_number_of_distinct_color_sets = [](auto& obj) { return obj.number_of_distinct_color_sets(); };
    auto call_sum_of_all_distinct_color_set_lengths = [](auto& obj) { return obj.sum_of_all_distinct_color_set_lengths(); };
    auto call_space_breakdown = [](auto& obj) { return obj.space_breakdown(); };


    cout << "Node length k: " << SBWT.get_k() << endl;
    cout << "Number of k-mers: " << SBWT.number_of_kmers() << endl;
    cout << "Number of subsets in the SBWT data structure: " << SBWT.number_of_subsets() << endl;
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
        sbwt::SeqIO::NullStream ns;
        int64_t bytes = SBWT.serialize(ns);
        cout << "SBWT: " << human_readable_bytes(bytes) << endl;
        cout << "==" << endl;

    }

    if(do_unitigs){
        write_log("Extracting unitigs (this could take a while)", LogLevel::MAJOR);

        string unitigs_file = get_temp_file_manager().create_filename("unitigs-",".fna");
        throwing_ofstream unitigs_out(unitigs_file);
        sbwt::SeqIO::NullStream null_stream;

        auto call_extract_unitigs = [&](auto& obj) {
            UnitigExtractor<decltype(obj)> UE;
            UE.extract_unitigs(dbg, obj, unitigs_out.stream, false, null_stream, null_stream);
        };

        std::visit(call_extract_unitigs, coloring);
        unitigs_out.close();

        LL unitig_count = 0;
        sbwt::SeqIO::Reader<> sr(unitigs_file);
        LL min_unitig_len = 1e18;
        LL max_unitig_len = 0;
        LL unitig_len_sum = 0;
        while(true){
            LL len = sr.get_next_read_to_buffer();
            if(len == 0) break;
            unitig_count++;
            min_unitig_len = min(min_unitig_len, len);
            max_unitig_len = max(max_unitig_len, len);
            unitig_len_sum += len;
        }

        cout << "Min unitig length: " << min_unitig_len << endl;
        cout << "Max unitig length: " << max_unitig_len << endl;
        cout << "Avg unitig length: " << (double)unitig_len_sum / unitig_count << endl;

    }

    return 0;

}
