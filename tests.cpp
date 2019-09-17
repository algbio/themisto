#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include "BD_BWT_index.hh"
#include "mark_kmers.hh"
#include "input_reading.hh"
#include "globals.hh"
#include "BOSS.hh"
#include "KallistoLite.hh"
#include <stdlib.h> 

using namespace std;

typedef int64_t LL; // long long
// The BWT will be the concatenation of all reads in sample A and then sample B with all reads separated by read_separator.

int main(){

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

    //test_build_metapangenome("./temp");

    test_pseudoalign("./temp");
    test_coloring();
    test_BOSS();

    

}