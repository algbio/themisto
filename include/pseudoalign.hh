#pragma once

#include <string>
#include "sbwt/SBWT.hh"
#include "new_coloring.hh"
#include "SeqIO.hh"
#include "variants.hh"

using namespace std;
using namespace sbwt;

typedef uint32_t color_t;

// Stores the intersection into buf1 and returns the number of elements in the
// intersection (does not resize buf1). Buffer elements must be sorted.
// Assumes all elements in a buffer are distinct
LL intersect_buffers(vector<color_t>& buf1, LL buf1_len, vector<color_t>& buf2, LL buf2_len){

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

// Stores the union into result_buf and returns the number of elements in the
// union (does not resize result_buf). Buffers elements must be sorted.
// Assumes all elements in a buffer are distinct. result_buf must have enough
// space to accommodate the union
LL union_buffers(vector<color_t>& buf1, LL buf1_len, vector<color_t>& buf2, LL buf2_len, vector<color_t>& result_buf){

    auto end = std::set_union(
                    buf1.begin(), buf1.begin() + buf1_len, 
                    buf2.begin(), buf2.begin() + buf2_len, 
                    result_buf.begin()
    );
    return end - result_buf.begin();
}

void pseudoalign(const plain_matrix_sbwt_t& SBWT, const Coloring& coloring, int64_t n_threads, std::string inputfile, std::string outputfile, bool reverse_complements, int64_t buffer_size, bool gzipped, bool sorted_output){

    SeqIO::Reader<> reader(inputfile);
    while(true) { 
        LL len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        vector<int64_t> colex_ranks = SBWT.streaming_search(reader.read_buf, len);
        Color_Set intersection;
        bool first_nonempty_found = false;
        for(LL colex : colex_ranks){
            if(colex >= 0){ // k-mer found Found
                Color_Set cs = coloring.get_color_set(colex);
                if(cs.size() > 0){
                    if(!first_nonempty_found){
                        intersection = cs;
                        first_nonempty_found = true;
                    } else{
                        intersection = intersection.intersection(cs);
                    }
                }
            }
        }
        // TODO: Print result
    }

}