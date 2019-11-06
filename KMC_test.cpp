#include "stdafx.h"
#include "external_memory_sort.hh"
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <functional>
#include "timer.h"
#include "globals.hh"
#include "BOSS.hh"

using namespace std;

// list_kmers is defined in KMC_wrapper.cpp
extern void list_kmers(int64_t k, int64_t ram_gigas, int64_t n_threads, 
                       string fastafile, string outfile, string tempdir);


struct Kmer{

    public:

    static const int MAX_SIZE = 256;

    char kmer[MAX_SIZE+1];

    void build_from_string(std::string S){
        assert(S.size() <= MAX_SIZE);
        for(int i = 0; i < S.size(); i++) kmer[i] = S[i];
        kmer[S.size()] = 0;
    }

    std::string to_string(){
        return kmer;
    }
};

struct colex_sort{
    bool operator()(const Kmer &a, const Kmer &b) const{
        std::string S1 = a.kmer;
        std::string S2 = b.kmer;
        std::reverse(S1.begin(), S1.end());
        std::reverse(S2.begin(), S2.end());
        return S1 < S2;
    }

    static Kmer min_value(){
        Kmer empty;
        empty.kmer[0] = 0;
        return empty;
    }

    static Kmer max_value(){
        Kmer maxval;
        maxval.kmer[0] = '~'; // ASCII code 126
        maxval.kmer[1] = 0;
        return maxval;
    }
};

// Adds cyclic kmers that cross some separator.
// APPENDS to outfile, does not overwrite
void add_extra_kmers(string fastafile, string outfile, LL k){
    throwing_ofstream kmerfile(outfile, std::ofstream::out | std::ofstream::app);
    FASTA_reader fr(fastafile);
    string concat;
    while(!fr.done()){
        string read = fr.get_next_query_stream().get_all();
        if(read.size() < 2*k){
            concat += read_separator + read;
        } else{
            concat += read_separator + read.substr(0,k) + read.substr(read.size()-k,k);
        }
    }
    concat += BD_BWT_index<>::END;
    for(string kmer : get_all_distinct_cyclic_kmers(concat,k)){
        bool good = false;
        for(char c : kmer) if(c == read_separator || c == BD_BWT_index<>::END){
            good = true;
        }
        if(good) kmerfile << kmer << "\n";
    }
}

void list_and_sort(string fastafile, string outfile, LL k){
    LL ram_gigas = 1;

    // Call KMC
    string kmc_outfile = temp_file_manager.get_temp_file_name("kmc_out");
    list_kmers(k, ram_gigas, 2, fastafile, kmc_outfile, temp_file_manager.get_dir());

    // Add kmers crossing separators and cyclic k-mers
    add_extra_kmers(fastafile, kmc_outfile, k);

    // Sort
    colex_sort cmp;
    external_memory_sort<Kmer, colex_sort>(kmc_outfile, outfile, cmp, 1024*1024*2014*ram_gigas);
}


void test_list_kmers(){
    Boss_Tester bt;

    // Generate testcases
    vector<pair<vector<string >, LL> > testcases; // pairs: (reads, k)
    LL n_reads = 10;
    for(LL read_length = 1; read_length <= 128; read_length *= 2){
        for(LL k = 1; k <= 256; k *= 2){
            vector<string> reads;
            for(LL i = 0; i < n_reads; i++){
                reads.push_back(get_random_dna_string(read_length));
            }
            testcases.push_back({reads,k});
        }
    }

    for(auto tcase : testcases){
        vector<string> reads = tcase.first;
        LL k = tcase.second;
        set<string> kmer_set = get_all_distinct_cyclic_kmers(reads,k);
        vector<string> colex_kmers(kmer_set.begin(), kmer_set.end());
        sort(colex_kmers.begin(), colex_kmers.end(), colex_compare);

        // Write reads out in fasta
        string fastafile = temp_file_manager.get_temp_file_name("fasta");
        throwing_ofstream fasta_out(fastafile);
        
        for(string read : reads){
            fasta_out << ">\n" << read << "\n";
        }
        fasta_out.close();

        // Call kmer listing
        string sorted_out = temp_file_manager.get_temp_file_name("sorted_kmers");
        list_and_sort(fastafile, sorted_out, k);

        // Read result
        vector<string> listed_kmers;
        string line;
        throwing_ifstream listed_in(sorted_out);
        while(listed_in.getline(line)){
            listed_kmers.push_back(line);
        }

        // Check
        assert(listed_kmers.size() == colex_kmers.size());
        for(LL i = 0; i < colex_kmers.size(); i++){
            cout << listed_kmers[i] << " " << colex_kmers[i] << endl;
            assert(listed_kmers[i] == colex_kmers[i]);
        }


    }

}


int main2(int argc, char** argv){

    set_temp_dir("temp");
    test_list_kmers();
    
}

int main(int argc, char** argv){
    try{
        return main2(argc, argv);
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    }
}