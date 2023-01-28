#include <iostream>
#include <unordered_map>
#include <string>
#include <set>
#include <exception>
#include <memory>
#include "SeqIO/SeqIO.hh"
#include "Kmer.hh"
#include "cxxopts.hpp"

typedef Kmer<32> kmer_t;

using namespace std;

vector<string> readlines(string filename){
    vector<string> lines;
    string line;
    throwing_ifstream in(filename);
    while(getline(in.stream,line)){
        lines.push_back(line);
    }
    return lines;
}

vector<int64_t> read_colorfile(string filename){
    vector<string> lines = readlines(filename);
    vector<int64_t> colors;
    for(const string& line : lines){
        colors.push_back(stoll(line));
    }
    return colors;
}

// Assumes gzipped data
vector<string> read_sequences(vector<string> filenames){
    vector<string> seqs;
    for(const string& f : filenames){
        SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in(f);
        while(true){
            int64_t len = in.get_next_read_to_buffer();
            if(len == 0) break;
            seqs.push_back(string(in.read_buf));
        }
    }
    return seqs;
}

// Returns true iff S[start .. start + k - 1] has only characters A,C,G and T
bool is_good_kmer(const char* S, int64_t start, int64_t k){
    for(int64_t i = start; i < start + k; i++){
        if(S[i] != 'A' && S[i] != 'C' && S[i] != 'G' && S[i] != 'T') return false;
    }
    return true;
}

unique_ptr<unordered_map<kmer_t, set<int64_t>>> build_kmer_map(const vector<string>& seqs, const vector<int64_t> colors, int64_t k){
    if(seqs.size() != colors.size()){
        cerr << seqs.size() << " " << colors.size() << endl;
        throw std::runtime_error("seqs.size() != colors.size()");
    }
    auto kmer_to_color_set = std::make_unique<unordered_map<kmer_t, set<int64_t>>>();
    for(int64_t seq_idx = 0; seq_idx < seqs.size(); seq_idx++){
        string S = seqs[seq_idx];
        for(char& c : S) c = toupper(c);

        int64_t color = colors[seq_idx];
        for(int64_t i = 0; i < (int64_t)S.size()-k+1; i++){
            if(is_good_kmer(S.c_str(), i, k)){
                kmer_t x(S.c_str() + i, k);
                (*kmer_to_color_set)[x].insert(color);
            }
        }
    }
    return kmer_to_color_set;
}

int main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Themisto reference implementation for testing");

    options.add_options()
        ("k", "The k of the k-mers.", cxxopts::value<int64_t>())
        ("i", "A text file with list of input files, one per line (.fna.gz files).", cxxopts::value<string>())
        ("c", "A file containing one integer color per sequence, one color per line.", cxxopts::value<string>())
        ("o", "The output directory..", cxxopts::value<string>())
    ;

    if(argc == 1){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    auto opts = options.parse(argc, argv);
    int64_t k = opts["k"].as<int64_t>();
    string in_file_list = opts["i"].as<string>();
    string color_file = opts["c"].as<string>();
    string out_dir = opts["o"].as<string>();

    vector<string> seqs = read_sequences(readlines(in_file_list));
    vector<int64_t> colors = read_colorfile(color_file);

    cerr << "Indexing..." << endl;
    unique_ptr<unordered_map<kmer_t, set<int64_t>>> kmer_to_color_set_ptr = build_kmer_map(seqs, colors, k);
    cerr << "Indexing done" << endl;
}