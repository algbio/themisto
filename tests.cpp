#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include "input_reading.hh"
#include "globals.hh"
#include "BOSS.hh"
#include "Themisto.hh"
#include <stdlib.h> 
#include "EM_sort_tests.hh"
#include "BOSS_tests.hh"
#include "version.h"
#include "Kmer_tests.hh"

using namespace std;

typedef int64_t LL; // long long

int main2(int argc, char** argv){
    // If no CLI parameters are given, run all tests. Otherwise run
    // only tests that are specified
    bool run_pseudoalign = argc == 1;
    bool run_BOSS = argc == 1;
    bool run_coloring = argc == 1;
    bool run_EM_sort = argc == 1;
    bool run_kmer_tests = argc == 1;

    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "-k") run_kmer_tests = true;
        if(string(argv[i]) == "-p") run_pseudoalign = true;
        if(string(argv[i]) == "-b") run_BOSS = true;
        if(string(argv[i]) == "-c") run_coloring = true;
        if(string(argv[i]) == "-s") run_EM_sort = true;
    }

    if(system("mkdir -p temp") != 0){
        cerr << "Error creating directory ./temp" << endl;
    }

    if(system("mkdir -p test_data") != 0){
        cerr << "Error creating directory ./test_data" << endl;
    }

    if(system("mkdir -p test_out") != 0){
        cerr << "Error creating directory ./test_out" << endl;
    }

    temp_file_manager.set_dir("temp");

    disable_logging();

    if(run_kmer_tests) kmer_tests();
    if(run_pseudoalign) test_pseudoalign("./temp");
    if(run_BOSS) test_BOSS();
    if(run_coloring) test_coloring();
    if(run_EM_sort) test_EM_sort();    
    write_log("All tests passed");
    
    return 0;

}

int main(int argc, char** argv){
    write_log("Themisto-" + std::string(THEMISTO_BUILD_VERSION));
    write_log("Maximum k-mer length: " + std::to_string(KMER_MAX_LENGTH));
    if(KMER_MAX_LENGTH != 255){
        write_log("Error: tests must be compiled with KMER_MAX_LENGTH=255");
        return 1;
    }
    (void)argc; (void)argv; // not used
    try{
        return main2(argc, argv);
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    }
}
