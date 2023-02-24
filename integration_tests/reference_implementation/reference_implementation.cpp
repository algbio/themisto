#include <iostream>
#include <unordered_map>
#include <string>
#include <set>
#include <exception>
#include <memory>
#include "SeqIO/SeqIO.hh"
#include "Kmer.hh"
#include "cxxopts.hpp"
#include "utils.hh"
#include "KmerIndex.hh"

typedef Kmer<32> kmer_t;
typedef Kmer<128> big_kmer_t;

using namespace std;

template<typename kmer_t> 
void run_queries(const vector<string>& seqs, const vector<int64_t> colors, const vector<string>& queries, int64_t k, double threshold, bool ignore_unknown, bool revcomps, string out_file){
    cerr << "Indexing..." << endl;
    KmerIndex<kmer_t> index(seqs, colors, k);
    cerr << "Indexing done" << endl;

    Buffered_ofstream<> out(out_file);
    cerr << "Querying" << endl;
    index.threshold_pseudoalign(queries, out, threshold, ignore_unknown, revcomps);
    cerr << "Querying done" << endl;    
}

template<typename kmer_t> 
void run_dump_color_matrix(const vector<string>& seqs, const vector<int64_t> colors, int64_t k, string out_file){
    cerr << "Indexing..." << endl;
    KmerIndex<kmer_t> index(seqs, colors, k);
    cerr << "Indexing done" << endl;
    cerr << "Dumping the color matrix" << endl;
    Buffered_ofstream<> out(out_file);
    index.dump_color_matrix(out);
    cerr << "Done" << endl;
}
  

int query_main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Themisto reference implementation for testing");

    options.add_options()
        ("k", "The k of the k-mers.", cxxopts::value<int64_t>())
        ("i", "A text file with list of input files, one per line (.fna.gz files).", cxxopts::value<string>())
        ("q", "A file with the query sequences (.fna.gz).", cxxopts::value<string>())
        ("threshold", "The pseudoalignment threshold.", cxxopts::value<double>())
        ("include-unknown-kmers", "Whether to include unknown k-mers.", cxxopts::value<bool>()->default_value("false"))
        ("rc", "Whether to also check against reverse complements.", cxxopts::value<bool>()->default_value("false"))
        ("o", "The output filename.", cxxopts::value<string>())
    ;

    if(argc == 1){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    auto opts = options.parse(argc, argv);
    int64_t k = opts["k"].as<int64_t>();
    string in_file_list = opts["i"].as<string>();
    string query_file = opts["q"].as<string>();
    string out_file = opts["o"].as<string>();
    double threshold = opts["threshold"].as<double>();
    double ignore_unknown =! opts["include-unknown-kmers"].as<bool>();
    bool revcomps = opts["rc"].as<bool>();
    
    vector<string> seqs = read_sequences(readlines(in_file_list));
    vector<string> queries = read_sequences({query_file});
    vector<int64_t> colors;
    for(int64_t i = 0; i < seqs.size(); i++) 
        colors.push_back(i); // Sequence colors

    if(k <= 32){
        run_queries<Kmer<32>>(seqs, colors, queries, k, threshold, ignore_unknown, revcomps, out_file);
    } else if(k <= 128){
        run_queries<Kmer<128>>(seqs, colors, queries, k, threshold, ignore_unknown, revcomps, out_file);
    } else{
        throw std::runtime_error("k too large");
        return 1;
    }

    return 0;
}

int dump_color_martrix_main(int argc, char** argv){
    cxxopts::Options options(argv[0], "Themisto reference implementation for testing");

    options.add_options()
        ("k", "The k of the k-mers.", cxxopts::value<int64_t>())
        ("i", "A text file with list of input files, one per line (.fna.gz files).", cxxopts::value<string>())
        ("rc", "Whether to add reverse complemets to the index.", cxxopts::value<bool>()->default_value("false"))
        ("file-colors", "", cxxopts::value<bool>()->default_value("false"))
        ("sequence-colors", "", cxxopts::value<bool>()->default_value("false"))
        ("manual-colors", "A manual colorfile, one per sequence", cxxopts::value<string>()->default_value(""))
        ("o", "The output filename.", cxxopts::value<string>())
    ;

    if(argc == 1){
        std::cerr << options.help() << std::endl;
        return 1;
    }

    auto opts = options.parse(argc, argv);
    int64_t k = opts["k"].as<int64_t>();
    string in_file_list = opts["i"].as<string>();
    string out_file = opts["o"].as<string>();
    bool revcomps = opts["rc"].as<bool>();
    bool file_colors = opts["file-colors"].as<bool>();
    bool sequence_colors = opts["sequence-colors"].as<bool>();
    string colorfile = opts["manual-colors"].as<string>();
    bool manual_colors = colorfile != "";

    vector<string> in_files = readlines(in_file_list);
    vector<string> seqs = read_sequences(in_files);

    vector<int64_t> colors;
    if(file_colors){
        for(int64_t color = 0; color < in_files.size(); color++){
            int64_t n_seqs = count_sequences(in_files[color]);
            for(int64_t i = 0; i < n_seqs; i++) colors.push_back(color);
        }
    } else if(sequence_colors){
        for(int64_t i = 0; i < seqs.size(); i++) colors.push_back(i);
    } else if(manual_colors){
        colors = read_colorfile(colorfile);
    } else{
        cerr << "Error: coloring mode not specified" << endl;
        return 1;
    }

    if(revcomps){
        int64_t n = seqs.size(); // Need to save this size to a local variable or else the loop below goes on forever
        for(int64_t i = 0; i < n; i++){
            seqs.push_back(get_rc(seqs[i]));
            colors.push_back(colors[i]);
        }
    }

    if(k <= 32){
        run_dump_color_matrix<Kmer<32>>(seqs, colors, k, out_file);
    } else if(k <= 128){
        run_dump_color_matrix<Kmer<128>>(seqs, colors, k, out_file);
    } else{
        throw std::runtime_error("k too large");
        return 1;
    }
    return 0;
}

int main(int argc, char** argv){
    if(argc == 1){
        cerr << "Commands: query, dump-color-matrix" << endl;
        return 1;
    }
    if(argv[1] == string("query"))
        return query_main(argc-1, argv+1);
    else if(argv[1] == string("dump-color-matrix"))
        return dump_color_martrix_main(argc-1, argv+1);
    else{
        cerr << "Unkown command " << argv[1] << endl;
    }
}
