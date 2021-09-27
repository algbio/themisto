#include "Themisto.hh"
#include "input_reading.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"

using namespace std;

char get_rc(char c){
    switch(c){
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return c;
    }
}

vector<string> read_lines(string filename){
    check_readable(filename);
    vector<string> lines;
    throwing_ifstream in(filename);
    string line;
    while(in.getline(line)){
        lines.push_back(line);
    }
    return lines;
}

struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const
    {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

int main2(int argc, char** argv){

    // Legacy support: transform old option format --outfile --out-file
    string legacy_support_fix = "--out-file";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--outfile") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "");

    options.add_options()
        ("i,index-dir", "Directory where the index will be built. Always required, directory must exist before running.", cxxopts::value<string>())
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    // todo: positional argument

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    string index_dir = opts["index-dir"].as<string>();
    check_true(index_dir != "", "Index directory not set");
    check_dir_exists(index_dir);

    string outfile = opts["out-file"].as<string>();
    check_writable(outfile);

    write_log("Starting");

    //temp_file_manager.set_dir(""); // There should not be any temp files

    write_log("Loading the index");    
    Themisto themisto;
    themisto.load_boss(index_dir + "/boss-");
    themisto.load_colors(index_dir + "/coloring-");

    write_log("Computing the counts");
    throwing_ofstream out(outfile);

    unordered_map<pair<LL, LL>, LL, hash_pair> pair_counts;
    vector<LL> color_counts(themisto.coloring.n_colors);

    BOSS<sdsl::bit_vector>& boss = themisto.boss;
    Progress_printer pp(themisto.boss.number_of_nodes(), 100);
    for(LL v = 0; v < themisto.boss.number_of_nodes(); v++){
        set<LL> colorset = themisto.coloring.get_colorset(v, themisto.boss);
        vector<LL> colorvec(colorset.begin(), colorset.end());
        for(LL i = 0; i < colorvec.size(); i++){
            color_counts[colorvec[i]]++;
            for(LL j = i+1; j < colorvec.size(); j++){
                LL c1 = colorvec[i];
                LL c2 = colorvec[j];
                if(c1 > c2) std::swap(c1,c2);
                pair_counts[{c1,c2}]++;
            }
        }
        pp.job_done();
    }

    write_log("Writing output");
    for(auto keyval : pair_counts){
        LL c1 = keyval.first.first; // Color 1
        LL c2 = keyval.first.second; // Color 2
        LL n_c1_c2 = keyval.second;
        LL n_c1 = color_counts[c1];
        LL n_c2 = color_counts[c2];
        out << c1 << " " << c2 << " " << n_c1 << " " << n_c2 << " " << n_c1_c2 <<  "\n";
    }

    write_log("Finished");

    return 0;
}

int main(int argc, char** argv){
    write_log("Themisto-" + std::string(THEMISTO_BUILD_VERSION));
    write_log("Maximum k-mer length (size of the de Bruijn graph node labels): " + std::to_string(KMER_MAX_LENGTH-1));
    try{
        return main2(argc, argv);
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    }
}
