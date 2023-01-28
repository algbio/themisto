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

using namespace std;

int main(int argc, char** argv){

    cxxopts::Options options(argv[0], "Themisto reference implementation for testing");

    options.add_options()
        ("k", "The k of the k-mers.", cxxopts::value<int64_t>())
        ("i", "A text file with list of input files, one per line (.fna.gz files).", cxxopts::value<string>())
        ("q", "A file with the query sequences (.fna.gz).", cxxopts::value<string>())
        ("threshold", "The pseudoalignment threshold.", cxxopts::value<double>())
        ("ignore-unknown-kmers", "Whether to ignore unknown k-mers.", cxxopts::value<bool>()->default_value("false"))
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
    double ignore_unknown = opts["ignore-unknown-kmers"].as<bool>();
    bool revcomps = opts["rc"].as<bool>();
    
    vector<string> seqs = read_sequences(readlines(in_file_list));
    vector<string> queries = read_sequences({query_file});
    vector<int64_t> colors;
    for(int64_t i = 0; i < seqs.size(); i++) 
        colors.push_back(i); // Sequence colors

    cerr << "Indexing..." << endl;
    KmerIndex index(seqs, colors, k);
    cerr << "Indexing done" << endl;

    Buffered_ofstream<> out(out_file);
    cerr << "Aligning" << endl;
    index.threshold_pseudoalign(queries, out, threshold, ignore_unknown, revcomps);
    cerr << "Aligning done" << endl;
}