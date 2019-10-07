#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <functional>
#include "timer.h"
#include "globals.hh"
#include "EM_sort.hh"
#include <cassert>

// Todo: remove this file from repository?

using namespace std;

// KMC_wrapper is defined in KMC_wrapper.cpp
extern void KMC_wrapper(int64_t k, int64_t ram_gigas, int64_t n_threads, 
                       string fastafile, string outfile, string tempdir);

// Also deletes leading empty lines as a side effect
void EM_delete_duplicate_lines(string infile, string outfile){

    ifstream in(infile);
    ofstream out(outfile);

    string prev;
    string line;
    while(getline(in, line)){
        if(line != prev) out << line << "\n";
        prev = line;
    }
}


int main(int argc, char** argv){

    string fastafile = argv[1];
    string tempdir = argv[2];
    LL RAM_gigas = stoll(argv[3]);
    LL n_threads = stoll(argv[4]);
    LL k = 33;

    set_temp_dir(tempdir);

    fastafile = fix_alphabet(fastafile);

    string KMC_file = temp_file_manager.get_temp_file_name("");

    write_log("Calling KMC");
    KMC_wrapper(k, RAM_gigas, n_threads, fastafile, KMC_file, temp_file_manager.get_dir());

    write_log("Sorting KMC");
    string KMC_sorted = temp_file_manager.get_temp_file_name("");
    EM_sort(KMC_file, KMC_sorted, colex_compare_cstrings, RAM_gigas * 1e9, 4, n_threads, EM_LINES);
    temp_file_manager.delete_file(KMC_file);

    write_log("Listing by brute");
    string brute_file = temp_file_manager.get_temp_file_name("");
    ofstream brute_out(brute_file);
    FASTA_reader fr(fastafile);
    while(!fr.done()){
        string seq = fr.get_next_query_stream().get_all();
        for(LL i = 0; i < seq.size(); i++){
            LL end = i + k - 1;
            if(end < seq.size()) brute_out << seq.substr(i,k) << "\n";
        }
    }
    brute_out.close();

    write_log("Sorting brute");
    string brute_sorted = temp_file_manager.get_temp_file_name("");
    EM_sort(brute_file, brute_sorted, colex_compare_cstrings, RAM_gigas * 1e9, 4, n_threads, EM_LINES);
    temp_file_manager.delete_file(brute_file);

    write_log("Deleting duplicate lines");
    string brute_uniques = temp_file_manager.get_temp_file_name("");
    EM_delete_duplicate_lines(brute_sorted, brute_uniques);
    temp_file_manager.delete_file(brute_sorted);

    write_log("Comparing KMC and brute");

    string brute_line, kmc_line;
    ifstream brute_in(brute_uniques);
    ifstream KMC_in(KMC_sorted);
    while(KMC_in || brute_in){
        getline(KMC_in, kmc_line);
        getline(brute_in, brute_line);
        if(!KMC_in) assert(!brute_in);
        if(!brute_in) assert(!KMC_in);
        if(kmc_line != brute_line){
            cout << "mismatch " << kmc_line << " " << brute_line << endl;
        }
    }

}

