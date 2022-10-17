#pragma once

#include "sdsl/int_vector.hpp"
#include <vector>
#include <algorithm>

class Fixed_Width_Int_Color_Set{

private:

    sdsl::int_vector<> v;

    // Number of bits required to represent x
    int64_t bits_needed(uint64_t x){
        int64_t ans = 0;
        while(x > 0){
            x >>= 1;
            ans++;
        }
        return max(ans, (int64_t) 1); // Need at least 1 bit even for zero
    }

    // Stores the intersection into buf1 and returns the number of elements in the
    // intersection (does not resize buf1). Buffer elements must be sorted.
    // Assumes all elements in a buffer are distinct
    int64_t intersect_buffers(vector<int64_t>& buf1, LL buf1_len, vector<int64_t>& buf2, int64_t buf2_len) const{

        int64_t i = 0, j = 0, k = 0;
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
    int64_t union_buffers(vector<int64_t>& buf1, LL buf1_len, vector<int64_t>& buf2, LL buf2_len, vector<int64_t>& result_buf) const{

        auto end = std::set_union(
                        buf1.begin(), buf1.begin() + buf1_len,
                        buf2.begin(), buf2.begin() + buf2_len,
                        result_buf.begin()
        );
        return end - result_buf.begin();
    }

public:

    Fixed_Width_Int_Color_Set() {}
    Fixed_Width_Int_Color_Set(const vector<int64_t>& colors){
        if(colors.size() > 0){
            int64_t max = *std::max_element(colors.begin(), colors.end());
            v = sdsl::int_vector<>(colors.size(), 0, bits_needed(max));
        }
    }

    bool empty() const{
        return v.size() > 0;
    }

    int64_t size() const{
        return v.size();
    }

    int64_t size_in_bits() const{
        return sdsl::size_in_bytes(v) * 8;
    }

    // Takes linear time
    bool contains(int64_t color) const{
        for(int64_t x : v) if(color == x) return true;
        return false;
    }

    Fixed_Width_Int_Color_Set intersection(const Fixed_Width_Int_Color_Set& other) const{
        vector<int64_t> A_vec = this->get_colors_as_vector();
        vector<int64_t> B_vec = other.get_colors_as_vector();
        int64_t size = intersect_buffers(A_vec, A_vec.size(), B_vec, B_vec.size());
        A_vec.resize(size);
        return Fixed_Width_Int_Color_Set(A_vec);
    }

    Fixed_Width_Int_Color_Set do_union(const Fixed_Width_Int_Color_Set& other) const{
        vector<int64_t> A_vec = this->get_colors_as_vector();
        vector<int64_t> B_vec = other.get_colors_as_vector();
        vector<int64_t> AB_vec(A_vec.size() + B_vec.size()); // Output buffer
        int64_t size = union_buffers(A_vec, A_vec.size(), B_vec, B_vec.size(), AB_vec);
        AB_vec.resize(size);
        return AB_vec;
    }

    vector<int64_t> get_colors_as_vector() const{
        return vector<int64_t>(v.begin(), v.end());
    }

    int64_t serialize(std::ostream& os){
        return v.serialize(os);
    }

    void load(std::istream& is){
        v.load(is);
    }
};
