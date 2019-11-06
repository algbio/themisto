#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include "BD_BWT_index.hh"
#include "input_reading.hh"
#include "globals.hh"
#include "BOSS.hh"
#include "KallistoLite.hh"
#include <stdlib.h> 
#include "EM_sort_tests.hh"

using namespace std;

typedef int64_t LL; // long long

int main2(int argc, char** argv){

    (void)argc; (void)argv; // not used

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

    //disable_logging();

    test_pseudoalign("./temp");
    test_BOSS();    
    test_coloring();
    test_EM_sort();    
    write_log("All tests passed");
    
    return 0;

}

int main(int argc, char** argv){
    (void)argc; (void)argv; // not used
    try{
        return main2(argc, argv);
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    }
}