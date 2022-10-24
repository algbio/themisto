#include "new_new_coloring.hh"

// See header for description
int64_t intersect_buffers(sdsl::int_vector<>& buf1, int64_t buf1_len, const sdsl::int_vector<>& buf2, int64_t buf2_start, int64_t buf2_len){

    int64_t i = 0, j = 0, k = 0;
    while(i < buf1_len && j < buf2_len){
        if(buf1[i] < buf2[buf2_start + j]) i++;
        else if(buf1[i] > buf2[buf2_start + j]) j++;
        else{
            buf1[k] = buf1[i];
            i++; j++; k++;
        }
    }
    return k;
}


// See header for description
int64_t union_buffers(vector<int64_t>& buf1, int64_t buf1_len, vector<int64_t>& buf2, int64_t buf2_len, vector<int64_t>& result_buf){

    auto end = std::set_union(
                    buf1.begin(), buf1.begin() + buf1_len,
                    buf2.begin(), buf2.begin() + buf2_len,
                    result_buf.begin()
    );
    return end - result_buf.begin();
}

// See header for description
int64_t bitmap_vs_bitmap_intersection(sdsl::bit_vector& A, int64_t A_size, const sdsl::bit_vector& B, int64_t B_start, int64_t B_size){
    int64_t n = min(A_size, B_size);
    int64_t words = n / 64;

    // Do 64-bit bitwise ands for the words that are in common
    for(int64_t w = 0; w < words; w++){
        A.set_int(w*64, A.get_int(w*64) & B.get_int(B_start + w*64));
    }

    // Do the rest one bit at a time. TODO: use a bit mask and do with just one and-operation
    for(int64_t i = words*64; i < n; i++){
        A[i] = A[i] && B[B_start + i];
    }

    return n;
}

// See header for description
int64_t array_vs_bitmap_intersection(sdsl::int_vector<>& iv, int64_t iv_size, const sdsl::bit_vector& bv, int64_t bv_start, int64_t bv_size){
    int64_t j = 0;
    for(int64_t i = 0; i < iv_size; i++){
        if(iv[i] >= bv_size) break;
        if(bv[bv_start + iv[i]]) iv[j++] = iv[i]; // Add to intersection modifying iv in-place
    }
    return j;
}

// See header for description
int64_t bitmap_vs_array_intersection(sdsl::bit_vector& bv, int64_t bv_size, const sdsl::int_vector<>& iv, int64_t iv_start, int64_t iv_size){
    int64_t iv_idx = 0;
    for(int64_t bv_idx = 0; bv_idx < bv_size; bv_idx++){
        int64_t iv_value = iv_idx < iv_size ? iv[iv_start + iv_idx] : -1;
        if(bv[bv_idx] == 1 && iv_value == bv_idx){
            bv[bv_idx] = 1; // Is in intersection
        } else{
            bv[bv_idx] = 0; // Is not in intersection
        }
        if(iv_value == bv_idx) iv_idx++;
    }
    return bv_size;
}

// See header for description
int64_t array_vs_array_intersection(sdsl::int_vector<>& A, int64_t A_len, const sdsl::int_vector<>& B, int64_t B_start, int64_t B_len){
    return intersect_buffers(A, A_len, B, B_start, B_len);
}

Color_Set_View::Color_Set_View(const Color_Set& cs) : start(cs.start), length(cs.length) {
        auto set_data_ptr = [&](auto& ptr){this->data_ptr = ptr;};
        std::visit(set_data_ptr, cs.data_ptr);
}