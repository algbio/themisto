#pragma once

#include "sdsl/bit_vectors.hpp"

class Bitmap_Or_Deltas_ColorSet{

// 4096 is a sampling parameter for random access. We don't really need random access so it is set to a high number.
typedef sdsl::enc_vector<coder::elias_delta, 4096> element_array_t;

private:

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

    Bitmap_Or_Deltas_ColorSet(const sdsl::bit_vector& bits) : is_bitmap(true) {
        bitmap = bits;
    }

    Bitmap_Or_Deltas_ColorSet(const element_array_t& elements) : is_bitmap(false) {
        element_array = elements;
    }


public:

    bool is_bitmap; // Is encoded as a bitmap or with gap encoding?

    sdsl::bit_vector bitmap;

    element_array_t element_array; // Todo: possibility of encoding deltas between non-existent colors

    Bitmap_Or_Deltas_ColorSet() : is_bitmap(true){}

    // Delete these constructors because they would be implicitly converted into an
    // element_array_t and then that constructor would be called, which we don't want
    Bitmap_Or_Deltas_ColorSet(const vector<std::uint64_t>& colors) = delete;
    Bitmap_Or_Deltas_ColorSet(const vector<std::int32_t>& colors) = delete;
    Bitmap_Or_Deltas_ColorSet(const vector<std::uint32_t>& colors) = delete;
    Bitmap_Or_Deltas_ColorSet(const vector<std::int16_t>& colors) = delete;
    Bitmap_Or_Deltas_ColorSet(const vector<std::uint16_t>& colors) = delete;
    Bitmap_Or_Deltas_ColorSet(const vector<std::int8_t>& colors) = delete;
    Bitmap_Or_Deltas_ColorSet(const vector<std::uint8_t>& colors) = delete;

    Bitmap_Or_Deltas_ColorSet(const vector<std::int64_t>& colors) {
        int64_t max_color = *std::max_element(colors.begin(), colors.end());
        
        element_array_t size_test(colors);
        if(sdsl::size_in_bytes(size_test)*8 > max_color+1){
            // Bitmap is smaller
            is_bitmap = true;
            sdsl::bit_vector bv(max_color+1, 0);
            for(int64_t x : colors) bv[x] = 1;
            bitmap = bv;
        } else{
            // Delta array is smaller
            is_bitmap = false;
            element_array = size_test;
        }
    }

    std::vector<int64_t> get_colors_as_vector() const {
        std::vector<int64_t> vec;
        if(is_bitmap){    
            for(int64_t i = 0; i < bitmap.size(); i++){
                if(bitmap[i]) vec.push_back(i);
            }
        } else{
            for(int64_t x : element_array){
                vec.push_back(x);
            }
        }
        return vec;
    }

    std::size_t size() const {
        if(is_bitmap){
            int64_t count = 0;
            for(bool b : bitmap) count += b;
            return count;
        }
        else return element_array.size();
    }


    std::size_t size_in_bits() const {
        return (sizeof(is_bitmap) + sdsl::size_in_bytes(bitmap) + sdsl::size_in_bytes(element_array)) * 8;
    }

    // This is O(1) time for dense color sets but O(set size) for sparse sets.
    bool contains(const int64_t color) const {
        if(color < 0) throw std::runtime_error("Called Color Set contains-method with a negative color id");
        if(is_bitmap) return color < bitmap.size() && bitmap[color];
        else{
            for(int64_t x : element_array) if(x == color) return true;
            return false;
        }
    }

    sdsl::bit_vector bitmap_vs_bitmap_intersection(const sdsl::bit_vector& A, const sdsl::bit_vector& B) const{
        int64_t n = min(A.size(), B.size());
        sdsl::bit_vector result(n, 0);
        int64_t words = n / 64;

        // Do 64-bit bitwise ands
        for(int64_t w = 0; w < words; w++){
            result.set_int(w*64, A.get_int(w*64) & B.get_int(w*64));
        }

        // Do the rest one bit at a time
        for(int64_t i = words*64; i < n; i++){
            result[i] = A[i] && B[i];
        }

        return result;    
    }

    sdsl::bit_vector bitmap_vs_bitmap_union(const sdsl::bit_vector& A, const sdsl::bit_vector& B) const{
        int64_t n = max(A.size(), B.size()); // Length of union
        sdsl::bit_vector result(n, 0);

        int64_t words = min(A.size(), B.size()) / 64; // Number of 64-bit words common to A and B

        // Do 64-bit bitwise ors
        for(int64_t w = 0; w < words; w++){
            result.set_int(w*64, A.get_int(w*64) | B.get_int(w*64));
        }

        // Do the rest one bit at a time
        for(int64_t i = words*64; i < n; i++){
            result[i] = (i < A.size() && A[i]) || (i < B.size() && B[i]);
        }

        return result;
    }

    element_array_t bitmap_vs_element_array_intersection(const sdsl::bit_vector& bm, const element_array_t& ea) const{
        vector<int64_t> new_elements;
        for(int64_t x : ea){
            if(x >= bm.size()) break;
            if(bm[x] == 1) new_elements.push_back(x);
        }
        return element_array_t(new_elements);
    }

    sdsl::bit_vector bitmap_vs_element_array_union(const sdsl::bit_vector& bm, const element_array_t& ea) const{
        if(ea.size() == 0) return bm;

        // Decode the integers in the element array
        vector<int64_t> elements;
        for(int64_t x : ea) elements.push_back(x);

        // Allocate space for the union
        int64_t union_size = max(elements.back()+1, (int64_t)bm.size());
        sdsl::bit_vector result(union_size, 0);

        // Add the bit map
        for(int64_t i = 0; i < bm.size(); i++) result[i] = bm[i];

        // Add the elements
        for(int64_t x : elements) result[x] = 1;

        return result;
    }    

    element_array_t element_array_vs_element_array_intersection(const element_array_t& A, const element_array_t& B) const{
        vector<int64_t> A_vec(A.begin(), A.end());
        vector<int64_t> B_vec(B.begin(), B.end());
        int64_t size = intersect_buffers(A_vec, A_vec.size(), B_vec, B_vec.size());
        A_vec.resize(size);
        return element_array_t(A_vec);
    }

    element_array_t element_array_vs_element_array_union(const element_array_t& A, const element_array_t& B) const{
        vector<int64_t> A_vec(A.begin(), A.end());
        vector<int64_t> B_vec(B.begin(), B.end());
        vector<int64_t> AB_vec(A_vec.size() + B_vec.size()); // Output buffer
        int64_t size = union_buffers(A_vec, A_vec.size(), B_vec, B_vec.size(), AB_vec);
        AB_vec.resize(size);
        return element_array_t(AB_vec);
    }

    Bitmap_Or_Deltas_ColorSet intersection(const Bitmap_Or_Deltas_ColorSet& c) const {
        if(is_bitmap && c.is_bitmap){
            return Bitmap_Or_Deltas_ColorSet(bitmap_vs_bitmap_intersection(bitmap, c.bitmap));
        } else if(is_bitmap && !c.is_bitmap){
            return Bitmap_Or_Deltas_ColorSet(bitmap_vs_element_array_intersection(bitmap, c.element_array));
        } else if(!is_bitmap && c.is_bitmap){
            return Bitmap_Or_Deltas_ColorSet(bitmap_vs_element_array_intersection(c.bitmap, element_array));
        } else{ // Element array vs element array
            return Bitmap_Or_Deltas_ColorSet(element_array_vs_element_array_intersection(element_array, c.element_array));
        }
    }

    // union is a reserved word in C++ so this function is called do_union
    Bitmap_Or_Deltas_ColorSet do_union(const Bitmap_Or_Deltas_ColorSet& c) const {
        if(is_bitmap && c.is_bitmap){
            return Bitmap_Or_Deltas_ColorSet(bitmap_vs_bitmap_union(bitmap, c.bitmap));
        } else if(is_bitmap && !c.is_bitmap){
            return Bitmap_Or_Deltas_ColorSet(bitmap_vs_element_array_union(bitmap, c.element_array));
        } else if(!is_bitmap && c.is_bitmap){
            return Bitmap_Or_Deltas_ColorSet(bitmap_vs_element_array_union(c.bitmap, element_array));
        } else{ // Element array vs element array
            return Bitmap_Or_Deltas_ColorSet(element_array_vs_element_array_union(element_array, c.element_array));
        }
    }

    std::size_t serialize(std::ostream& os) const {
        int64_t n_bytes_written = 0;

        char flag = is_bitmap;
        os.write(&flag, 1);
        n_bytes_written  += 1;

        n_bytes_written += bitmap.serialize(os);
        n_bytes_written += element_array.serialize(os);

        return n_bytes_written;
    }

    void load(std::ifstream& is) {
        char flag;
        is.read(&flag, 1);
        is_bitmap = flag;
        bitmap.load(is);
        element_array.load(is);
    }

};


class SDSL_Hybrid_Set {
    sdsl::hyb_vector<> data;

public:

    SDSL_Hybrid_Set(){}

    SDSL_Hybrid_Set(const vector<std::int64_t>& colors) {
        int64_t max_color = *std::max_element(colors.begin(), colors.end());
        sdsl::bit_vector vec(max_color+1, 0);
        for(int64_t x : colors) vec[x] = 1;
        data = sdsl::hyb_vector<>(vec);
    }


    std::vector<std::uint64_t> get_colors_as_vector() const {
        std::vector<std::uint64_t> vec;
        for(int64_t i = 0; i < data.size(); i++){
            if(data[i]) vec.push_back(i);
        }
        return vec;
    }

    std::size_t size() const {
        return get_colors_as_vector().size(); // TODO SLOW
    }


    std::size_t size_in_bits() const {
        return sdsl::size_in_bytes(data) * 8;
    }

    bool contains(const std::uint64_t n) const {
        return data[n];
    }

    SDSL_Hybrid_Set intersection(const SDSL_Hybrid_Set& c) const {
        // TODO: SLOW
        vector<std::int64_t> result;
        int64_t n = min(data.size(), c.data.size());
        for(int64_t i = 0; i < n; i++){
            if(data[i] && c.data[i])
                result.push_back(i);
        }
        return SDSL_Hybrid_Set(result);
    }

    // union is a reserved word in C++ so this function is called do_union
    SDSL_Hybrid_Set do_union(const SDSL_Hybrid_Set& c) const {
        // TODO: SUPER SLOW
        set<std::int64_t> result;
        for(std::int64_t x : get_colors_as_vector()){
            result.insert(x);
        }
        for(std::int64_t x : c.get_colors_as_vector()){
            result.insert(x);
        }

        vector<std::int64_t> vec(result.begin(), result.end());
        return SDSL_Hybrid_Set(vec);
    }

    std::size_t serialize(std::ostream& os) const {
        return data.serialize(os);
    }

    std::size_t load(std::ifstream& is) {
        data.load(is);
        return sdsl::size_in_bytes(data);
    }

};