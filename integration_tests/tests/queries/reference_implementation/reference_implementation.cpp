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

// Table mapping ascii values of characters to their reverse complements,
// lower-case to lower case, upper-case to upper-case. Non-ACGT characters
// are mapped to themselves.
static constexpr unsigned char rc_table[256] =
{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
56, 57, 58, 59, 60, 61, 62, 63, 64, 84, 66, 71, 68, 69, 70, 67, 72, 73,
74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 65, 85, 86, 87, 88, 89, 90, 91,
92, 93, 94, 95, 96, 116, 98, 103, 100, 101, 102, 99, 104, 105, 106, 107,
108, 109, 110, 111, 112, 113, 114, 115, 97, 117, 118, 119, 120, 121, 122,
123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137,
138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152,
153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182,
183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212,
213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227,
228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242,
243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255};

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

string get_rc(const string& S){
    string revc(S.rbegin(), S.rend());
    for(char& c : revc) c = rc_table[c];
    return revc;
}

template<typename output_stream_t>
void write_number_in_ascii(output_stream_t& out, int64_t x){
    string ascii = to_string(x);
    out.write(ascii.c_str(), ascii.size()); 
}

template<typename output_stream_t>
void threshold_pseudoalign(const unordered_map<kmer_t, set<int64_t>>& index, int64_t k, vector<string>& queries, output_stream_t& out, double threshold, bool ignore_unknown, bool revcomps){

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

                string kmer = query.substr(i, k);
                string kmer_revc = get_rc(kmer);

                kmer_t x(kmer);
                kmer_t x_revc(kmer_revc);

                set<int64_t> colors;

                bool known = false;
                if(index.find(x) != index.end()){ // Forward k-mer found
                    known = true;
                    for(int64_t color : index.at(x)){
                        colors.insert(color);        
                    }
                }

                if(revcomps && index.find(x_revc) != index.end()){ // Reverse complement k-mer found
                    known = true;
                    for(int64_t color : index.at(x_revc)){
                        colors.insert(color);        
                    }
                } 

                for(int64_t color : colors) counts[color]++;

                if(!known && ignore_unknown) n_relevant_kmers--;
            }
        }

        // Report results
        write_number_in_ascii(out, query_id);
        if(n_relevant_kmers == 0) continue;
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
    unique_ptr<unordered_map<kmer_t, set<int64_t>>> kmer_to_color_set = build_kmer_map(seqs, colors, k);
    cerr << "Indexing done" << endl;

    Buffered_ofstream<> out(out_file);
    cerr << "Aligning" << endl;
    threshold_pseudoalign(*kmer_to_color_set, k, queries, out, threshold, ignore_unknown, revcomps);
    cerr << "Aligning done" << endl;
}