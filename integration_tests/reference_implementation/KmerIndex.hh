#pragma once

#include <unordered_map>
#include <set>
#include "Kmer.hh"
#include "utils.hh"

using namespace std;

template<typename kmer_t>
class KmerIndex{

public:

    unordered_map<kmer_t, set<int64_t>> mapping; // kmer -> color set
    int64_t k;

    KmerIndex(const vector<string>& seqs, const vector<int64_t> colors, int64_t k) : k(k){
        if(seqs.size() != colors.size()){
            cerr << seqs.size() << " " << colors.size() << endl;
            throw std::runtime_error("seqs.size() != colors.size()");
        }

        for(int64_t seq_idx = 0; seq_idx < seqs.size(); seq_idx++){
            string S = seqs[seq_idx];
            int64_t color = colors[seq_idx];

            for(int64_t i = 0; i < (int64_t)S.size()-k+1; i++){
                if(is_good_kmer(S.c_str(), i, k)){
                    kmer_t x(S.c_str() + i, k);
                    mapping[x].insert(color);
                }
            }
        }
    }

    template<typename output_stream_t>
    void threshold_pseudoalign(const vector<string>& queries, output_stream_t& out, double threshold, bool ignore_unknown, bool revcomps){
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
                    if(mapping.find(x) != mapping.end()){ // Forward k-mer found
                        known = true;
                        for(int64_t color : mapping.at(x)){
                            colors.insert(color);        
                        }
                    }

                    if(revcomps && mapping.find(x_revc) != mapping.end()){ // Reverse complement k-mer found
                        known = true;
                        for(int64_t color : mapping.at(x_revc)){
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
                if((double) count / n_relevant_kmers >= threshold){
                    out.write(&space, 1);
                    write_number_in_ascii(out, color);
                }
            }
            out.write(&newline, 1);
            query_id++;
        }
    }

    template<typename output_stream_t>
    void dump_color_matrix(output_stream_t& out){
        vector<pair<kmer_t, set<int64_t>>> matrix; // matrix[i] will be the list of non-zero entries in i-th row of the the color matrix
        for(auto [kmer, colors] : mapping){
            matrix.push_back({kmer, colors});
        }

        std::sort(matrix.begin(), matrix.end()); // Sort by k-mer

        // Write out
        char space = ' ';
        char newline = '\n';
        for(auto [kmer, colors] : matrix){
            // Write the k-mer
            string kmer_str = kmer.to_string();
            out.write(kmer_str.c_str(), kmer_str.size());

            // Write the colors
            for(int64_t color : colors){
                out.write(&space, 1);
                write_number_in_ascii(out, color);
            }
            out.write(&newline, 1);
        }
    }

};