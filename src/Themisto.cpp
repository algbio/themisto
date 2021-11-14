#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <set>
#include <algorithm>
#include <filesystem>
#include "input_reading.hh"
#include "test_tools.hh"
#include "globals.hh"
#include "libwheeler/BOSS.hh"
#include "libwheeler/BOSS_builder.hh"
#include "Coloring.hh"
#include "WorkDispatcher.hh"
#include "EM_sort.hh"
#include "KMC_wrapper.hh"
#include "Themisto.hh"

using namespace std;

// Stores the intersection into buf1 and returns the number of elements in the
// intersection (does not resize buf1). Buffer elements must be sorted.
// Assumes all elements in a buffer are distinct
LL intersect_buffers(vector<LL>& buf1, LL buf1_len, vector<LL>& buf2, LL buf2_len){

    LL i = 0, j = 0, k = 0;
    while(i < buf1_len && j < buf2_len){
        if(buf1[i] < buf2[j]) i++;
        else if(buf1[i] > buf2[j]) j++;
        else{
            buf1[k] = buf1[i];
            i++; j++; k++;
        }
    }
    return k;

}

// Stores the union into buf1 and returns the number of elements in the
// intersection (does not resize buf1). Buffers elements must be sorted.
// Assumes all elements in a buffer are distinct. result_buf must have enough
// space to accommodate the union
LL union_buffers(vector<LL>& buf1, LL buf1_len, vector<LL>& buf2, LL buf2_len, vector<LL>& result_buf){

    auto end = std::set_union(
                    buf1.begin(), buf1.begin() + buf1_len, 
                    buf2.begin(), buf2.begin() + buf2_len, 
                    result_buf.begin()
    );
    return end - result_buf.begin();
}
