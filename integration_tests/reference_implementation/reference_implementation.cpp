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

template<typename output_stream_t>
void write_number_in_ascii(output_stream_t& out, int64_t x){
    string ascii = to_string(x);
    out.write(ascii.c_str(), ascii.size()); 
}

template<typename output_stream_t>
void threshold_pseudoalign(const unordered_map<kmer_t, set<int64_t>>& index, int64_t k, vector<string>& queries, output_stream_t& out, double threshold, bool ignore_unknown){

    char space = ' ';
    char newline = '\n';
    int64_t query_id = 0;
    for(const string& query : queries){
        // Look up k-mers
        unordered_map<int64_t, int64_t> counts; // color -> count
        int64_t n_relevant_kmers = 0;
        for(int64_t i = 0; i < (int64_t)query.size()-k+1; i++){ // For all k-mers
            if(is_good_kmer(query.c_str(), i, k)){ // Must have only DNA-characters
                n_relevant_kmers++;
                kmer_t x(query.c_str() + i, k);
                if(index.find(x) != index.end()){ // Known k-mer
                    for(int64_t color : index.at(x)){
                        counts[color]++;
                    }
                } else{ // Unknown k-mer
                    if(ignore_unknown) n_relevant_kmers--;
                }
            }
        }

        // Report results
        write_number_in_ascii(out, query_id);
        for(auto [color, count] : counts){
            if((double) count / n_relevant_kmers > threshold){
                out.write(&space, 1);
                write_number_in_ascii(out, color);
            }
        }
        out.write(&newline, 1);
        query_id++;
    }
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
    unique_ptr<unordered_map<kmer_t, set<int64_t>>> kmer_to_color_set = build_kmer_map(seqs, colors, k);
    cerr << "Indexing done" << endl;

    {Buffered_ofstream<> out(out_dir + "/theshold-ignore-0.01.txt");
    threshold_pseudoalign(*kmer_to_color_set, k, seqs, out, 0.01, true);}
    {Buffered_ofstream<> out(out_dir + "/theshold-ignore-0.1.txt");
    threshold_pseudoalign(*kmer_to_color_set, k, seqs, out, 0.1, true);}
    {Buffered_ofstream<> out(out_dir + "/theshold-ignore-0.5.txt");
    threshold_pseudoalign(*kmer_to_color_set, k, seqs, out, 0.5, true);}
    {Buffered_ofstream<> out(out_dir + "/theshold-ignore-0.9.txt");
    threshold_pseudoalign(*kmer_to_color_set, k, seqs, out, 0.9, true);}       
    {Buffered_ofstream<> out(out_dir + "/theshold-ignore-0.1.01.txt");
    threshold_pseudoalign(*kmer_to_color_set, k, seqs, out, 1.01, true);}

    {Buffered_ofstream<> out(out_dir + "/theshold-0.01.txt");
    threshold_pseudoalign(*kmer_to_color_set, k, seqs, out, 0.01, false);}
    {Buffered_ofstream<> out(out_dir + "/theshold-0.1.txt");
    threshold_pseudoalign(*kmer_to_color_set, k, seqs, out, 0.1, false);}
    {Buffered_ofstream<> out(out_dir + "/theshold-0.5.txt");
    threshold_pseudoalign(*kmer_to_color_set, k, seqs, out, 0.5, false);}
    {Buffered_ofstream<> out(out_dir + "/theshold-0.9.txt");
    threshold_pseudoalign(*kmer_to_color_set, k, seqs, out, 0.9, false);}
    {Buffered_ofstream<> out(out_dir + "/theshold-0.1.01.txt");
    threshold_pseudoalign(*kmer_to_color_set, k, seqs, out, 1.01, false);}
}