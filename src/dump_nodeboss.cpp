#include "Themisto.hh"
#include "input_reading.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"

using namespace std;

map<char, vector<bool>> get_nodeboss_bit_vectors(BOSS<sdsl::bit_vector>& boss){
    // Mark the colex-smallest edge to each node.

    map<char, vector<bool>> BWTs;
    for(char c : boss.get_alphabet())
        BWTs[c].resize(boss.number_of_nodes());

    LL node_id = -1; // Wheeler rank
    LL edge_id = -1; // Wheeler rank
    LL edge_count_to_current_node = 0;
    Progress_printer pp(boss.indegs_size(), 100);

    for(LL i = 0; i < boss.indegs_size(); i++){
        if(boss.indegs_at(i) == 1){
            node_id++;
            edge_count_to_current_node = 0;
        }
        else{
            // Process edge
            edge_id++;
            edge_count_to_current_node++;
            char label = boss.incoming_character(node_id);

            if(edge_count_to_current_node == 1){
                LL source_node = boss.edge_source(edge_id);
                BWTs[label][source_node] = 1;
            }
        }
        pp.job_done();
    }
    return BWTs;
}

void dump_bit_vector(const vector<bool>& vec, string filename){
    throwing_ofstream out(filename, ios::binary);
    int64_t size = vec.size();
    out.stream.write((char*)&size, sizeof(size));

    // Pad to a multiple of 8
    vector<bool> padded = vec;
    while(padded.size() % 8 != 0) padded.push_back(0);
    int64_t n_bytes = padded.size() / 8;

    for(int64_t i = 0; i < n_bytes; i++){
        uint8_t byte = (padded[8*i + 0] << 7) | 
                       (padded[8*i + 1] << 6) | 
                       (padded[8*i + 2] << 5) | 
                       (padded[8*i + 3] << 4) |
                       (padded[8*i + 4] << 3) |
                       (padded[8*i + 5] << 2) |
                       (padded[8*i + 6] << 1) |
                       (padded[8*i + 7] << 0);
        out.stream.write((char*)&byte, 1);
    }
}

vector<LL> get_C_array(map<char, vector<bool>>& BWTs){
    vector<LL> C(4);
    LL char_index = 0;
    LL cumul_count = 1; // Starting at 1 because there is a ghost dollar to the root node.
    for(auto X : BWTs){
        C[char_index] = cumul_count;
        vector<bool>& vec = X.second;
        for(LL i = 0; i < vec.size(); i++)
            cumul_count += vec[i];
        char_index++;
    }
    return C;
}

void dump_C_array(const vector<LL>& C, string filename){
    throwing_ofstream out(filename, ios::binary);
    int64_t size = C.size();
    out.stream.write((char*)&size, sizeof(size));

    for(LL i = 0; i < C.size(); i++){
        out.stream.write((char*)&C[i], sizeof(LL));
    }
}

LL nodeboss_search(const map<char, vector<LL>>& rank1_answers, const vector<LL>& C, const string& kmer){
    LL node_left = 0;
    LL node_right = rank1_answers.at('A').size()-1; // number of nodes - 1
    for(LL i = 0; i < kmer.size(); i++){
        cerr << node_left << " " << node_right << endl;
        char c = kmer[i];
        node_left = C[c] + rank1_answers.at(c)[node_left];
        node_right = C[c] + rank1_answers.at(c)[node_right+1] - 1;
        if(node_left > node_right) return -1; // Not found
    }
    cerr << node_left << " " << node_right << endl;
    if(node_left != node_right){
        cerr << "node_left != node_right" << endl;
        exit(1);
    }
    return node_left;
}


void run_test(const BOSS<sdsl::bit_vector>& boss, const map<char, vector<bool>>& BWTs, const vector<LL>& C){

    // Precompute rank answers
    write_log("Precomputing rank answers", LogLevel::MAJOR);
    map<char, vector<LL>> rank1_answers;
    for(char c : boss.get_alphabet()){
        rank1_answers[c].resize(BWTs.at(c).size());
        LL answer = 0;
        for(LL i = 0; i < rank1_answers[c].size(); i++){
            rank1_answers[c][i] = answer;
            answer += BWTs.at(c)[i];
        }
        cout << c << " " << rank1_answers[c] << endl;
    }

    // Search for all nodes
    for(LL v = 0; v < boss.number_of_nodes(); v++){
        string kmer = boss.get_node_label(v);
        if(kmer.size() < boss.get_k()) continue;
        LL colex = nodeboss_search(rank1_answers, C, kmer);
        cout << colex << " " << v << endl;
        if(colex != v){
            cout << "WRONG ANSWER" << endl;
            exit(1);
        }
    }

}

void dump_nodeboss(BOSS<sdsl::bit_vector>& boss, string out_prefix){

    write_log("Building the BWT bit vectors", LogLevel::MAJOR);
    map<char, vector<bool>> BWTs = get_nodeboss_bit_vectors(boss);

    write_log("Dumping bit BWT vectors to disk", LogLevel::MAJOR);
    for(char c : boss.get_alphabet()){
        cout << c << " ";; // DEBUG
        for(bool b : BWTs[c]) cout << b; // DEBUG
        cout << endl; // DEBUG
        string filename = out_prefix + "BWT_" + c;
        dump_bit_vector(BWTs[c], filename);
    }

    write_log("Computing the C array", LogLevel::MAJOR);
    vector<LL> C_array = get_C_array(BWTs);
    cout << "C array " << C_array << endl; // DEBUG

    write_log("Dumping C array to disk", LogLevel::MAJOR);
    dump_C_array(C_array, out_prefix + "_C");

    run_test(boss, BWTs, C_array);
    
}

int main(int argc, char** argv){

    set_log_level(LogLevel::MINOR);

    // Legacy support: transform old option format --outfile --out-file
    string legacy_support_fix = "--out-file";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--outfile") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "This program prints a NodeBOSS representation of the index. Not production quality code.");

    options.add_options()
        ("o,out-prefix", "Output file prefix.", cxxopts::value<string>())
        ("i,dbg-file", "The .tdbg file of an index.", cxxopts::value<string>())
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        cerr << "Usage examples: " << argv[0] << " -i in_prefix -o out_prefix --temp-dir temp" << endl;
        exit(1);
    }

    string out_prefix = opts["out-prefix"].as<string>();
    string index_dbg_file = opts["dbg-file"].as<string>();
    string temp_dir = opts["temp-dir"].as<string>();

    create_directory_if_does_not_exist(temp_dir);

    write_log("Starting", LogLevel::MAJOR);

    get_temp_file_manager().set_dir(temp_dir);

    write_log("Loading the DBG", LogLevel::MAJOR);
    Themisto themisto;
    themisto.load_boss(index_dbg_file);
    
    dump_nodeboss(themisto.boss, out_prefix);

    return 0;
}
