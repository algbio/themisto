#pragma once

#include "sdsl/int_vector.hpp"
#include <vector>
#include <algorithm>
#include <bit>

class Fixed_Width_Int_Color_Set{

private:

    sdsl::int_vector<> v;

    // Number of bits required to represent x
    int64_t bits_needed(uint64_t x){
        return max((int64_t)std::bit_width(x), (int64_t)1); // Need at least 1 bit (for zero)
    }

    // Stores the intersection into buf1 and returns the number of elements in the
    // intersection (does not resize buf1). Buffer elements must be sorted.
    // Assumes all elements in a buffer are distinct
    int64_t intersect_buffers(sdsl::int_vector<>& buf1, int64_t buf1_len, const sdsl::int_vector<>& buf2, int64_t buf2_len) const{

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
    int64_t union_buffers(const sdsl::int_vector<>& buf1, int64_t buf1_len, const sdsl::int_vector<>& buf2, int64_t buf2_len, sdsl::int_vector<>& result_buf) const{

        auto end = std::set_union(
                        buf1.begin(), buf1.begin() + buf1_len,
                        buf2.begin(), buf2.begin() + buf2_len,
                        result_buf.begin()
        );
        return end - result_buf.begin();
    }

public:

    Fixed_Width_Int_Color_Set() {}
    Fixed_Width_Int_Color_Set(const sdsl::int_vector<>& v) : v(v){}
    Fixed_Width_Int_Color_Set(const vector<int64_t>& colors){
        if(colors.size() > 0){
            int64_t max = *std::max_element(colors.begin(), colors.end());
            v = sdsl::int_vector<>(colors.size(), 0, bits_needed(max));
            for(int64_t i = 0; i < colors.size(); i++) v[i] = colors[i];
        }
    }

    bool empty() const{
        return v.size() == 0;
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
        sdsl::int_vector this_copy = v;
        int64_t size = intersect_buffers(this_copy, this_copy.size(), other.v, other.v.size());
        this_copy.resize(size);
        return Fixed_Width_Int_Color_Set(this_copy);
    }

    Fixed_Width_Int_Color_Set do_union(const Fixed_Width_Int_Color_Set& other) const{
        sdsl::int_vector<> result(this->v.size() + other.v.size()); // Space for the union
        int64_t size = union_buffers(v, v.size(), other.v, other.v.size(), result);
        result.resize(size);
        return Fixed_Width_Int_Color_Set(result);
    }

    vector<int64_t> get_colors_as_vector() const{
        return vector<int64_t>(v.begin(), v.end());
    }

    void push_colors_to_vector(std::vector<int64_t>& vec) const{
        for(int64_t color : v) vec.push_back(color);
    }


    int64_t serialize(std::ostream& os) const{
        return v.serialize(os);
    }

    void load(std::istream& is){
        v.load(is);
    }
};
