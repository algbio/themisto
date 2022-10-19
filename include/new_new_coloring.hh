#pragma once

#include <vector>
#include "sdsl/bit_vectors.hpp"
#include "sdsl/int_vector.hpp"
#include <variant>

using namespace std;

// This class stores a sequence of n non-negative integers with total sum N in (n + N) + o(n + N)
// bits of space, and can answer in constant time queries for sums of the first i stored integers.
class Succinct_Prefix_Sums{

    private:
    sdsl::bit_vector v;
    sdsl::rank_support_v5<0> rs;
    sdsl::select_support_mcl<1> ss;

    vector<bool> temp_v; // Used only during construction

    void append_unary_number(int64_t x, vector<bool>& v){
        v.push_back(1);
        for(int64_t i = 0; i < x; i++) v.push_back(0);
    }

    public:

    void add(int64_t x){
        append_unary_number(x, temp_v);
    }

    void finish_building(){
        temp_v.push_back(1); // End sentinel

        // Copy temp_v into v
        v = sdsl::bit_vector(temp_v.size(), 0);
        for(int64_t i = 0; i < v.size(); i++) v[i] = temp_v[i];

        sdsl::util::init_support(rs, &v);
        sdsl::util::init_support(ss, &v);
    }

    int64_t serialize(ostream& os) const{
        int64_t bytes_written = 0;
        bytes_written += v.serialize(os);
        bytes_written += rs.serialize(os);
        bytes_written += ss.serialize(os);

        return bytes_written;

        // Do not serialize temp_v
    }

    void load(istream& is){
        v.load(is);
        rs.load(is, &v);
        ss.load(is, &v);

        // Do not load temp_v
    }

    // Returns the sum of v[0..i), where v is the array of stored integers.
    // The parameter i can be from 0 to n inclusive, where n is the number of stored integers.
    int64_t sum(int64_t i) const{
        return rs.rank(ss.select(i+1));
    }
};

class New_Hybrid_Color_Set{
    
public:

    int64_t start;
    int64_t length;

    bool owns_memory;
    std::variant<const sdsl::bit_vector*, const  sdsl::int_vector<>*> data_ptr;

    bool is_bitmap(){
        return std::holds_alternative<const sdsl::bit_vector*>(data_ptr);
    }

    New_Hybrid_Color_Set(int64_t start, int64_t length, std::variant<const sdsl::bit_vector*, const sdsl::int_vector<>*> data_ptr)
        : start(start), length(length), owns_memory(false), data_ptr(data_ptr) {}

};

template<typename color_set_t> class Color_Set_Storage;

template<>
class Color_Set_Storage<New_Hybrid_Color_Set>{

    private:

    sdsl::bit_vector bitmap_concat;
    Succinct_Prefix_Sums bitmap_sizes;

    sdsl::int_vector<> deltas_concat;
    Succinct_Prefix_Sums deltas_sizes;

    sdsl::bit_vector is_bitmap_marks;
    sdsl::rank_support_v5<> is_bitmap_marks_rs;

    // Dynamic-length vectors used during construction only
    vector<bool> temp_bitmap_concat;
    vector<int64_t> temp_deltas_concat;
    vector<bool> temp_is_bitmap_marks;

    // Number of bits required to represent x
    int64_t bits_needed(uint64_t x){
        return max((int64_t)std::bit_width(x), (int64_t)1); // Need at least 1 bit (for zero)
    }

    sdsl::bit_vector to_sdsl_bit_vector(const vector<bool>& v){
        sdsl::bit_vector bv(v.size());
        for(int64_t i = 0; i < v.size(); i++) bv[i] = v[i];
        return bv;
    }

    public:

    New_Hybrid_Color_Set get_color_set_by_id(int64_t id) const{
        if(is_bitmap_marks[id]){
            int64_t bitmap_idx = is_bitmap_marks_rs.rank(id); // This many bitmaps come before this bitmap
            int64_t start = bitmap_sizes.sum(bitmap_idx);
            int64_t end = bitmap_sizes.sum(bitmap_idx+1); // One past the end
            std::variant<const sdsl::bit_vector*, const sdsl::int_vector<>*> data_ptr = &bitmap_concat;
            return New_Hybrid_Color_Set(start, end-start, data_ptr);
        } else{
            int64_t bitmap_idx = id - is_bitmap_marks_rs.rank(id); // Rank-0. This many delta arrays come before this bitmap
            int64_t start = deltas_sizes.sum(bitmap_idx);
            int64_t end = deltas_sizes.sum(bitmap_idx+1); // One past the end
            std::variant<const sdsl::bit_vector*, const sdsl::int_vector<>*> data_ptr = &deltas_concat;
            return New_Hybrid_Color_Set(start, end-start, data_ptr);
        }
    }

    // Need to call prepare_for_queries() after all sets have been added
    // Set must be sorted
    void add_set(const vector<int64_t>& set){

        int64_t max_element = *std::max_element(set.begin(), set.end());
        if(log2(max_element) * set.size() > max_element){
            // Dense -> bitmap

            // Add is_bitmap_mark
            temp_is_bitmap_marks.push_back(1);

            // Create bitmap
            vector<bool> bitmap(set.size(), 0);
            for(int64_t x : set) bitmap[x] = 1;
            for(bool b : bitmap) temp_bitmap_concat.push_back(b);

            // Store bitmap length
            bitmap_sizes.add(max_element+1);
        } else{
            // Sparse -> delta array

            // Add is_bitmap_mark
            temp_is_bitmap_marks.push_back(0);

            if(set.size() > 0){
                for(int64_t i = 0; i < set.size(); i++){
                    if(i == 0) temp_deltas_concat.push_back(set[0]);
                    else temp_deltas_concat.push_back(set[i] - set[i-1]);
                }
                
            }

            // Store array length
            deltas_sizes.add(set.size());
        }

    }


    // Call this after done with add_set
    void prepare_for_queries(){

        int64_t max_delta = *std::max_element(temp_deltas_concat.begin(), temp_deltas_concat.end());
        deltas_concat = sdsl::int_vector<>(temp_deltas_concat.size(), 0, bits_needed(max_delta));
        for(int64_t i = 0; i < temp_deltas_concat.size(); i++)
            deltas_concat[i] = temp_deltas_concat[i];

        bitmap_concat = to_sdsl_bit_vector(temp_bitmap_concat);
        is_bitmap_marks = to_sdsl_bit_vector(temp_is_bitmap_marks);
        
        bitmap_sizes.finish_building();
        deltas_sizes.finish_building();

        sdsl::util::init_support(is_bitmap_marks_rs, &is_bitmap_marks);

        // Free memory
        temp_deltas_concat.clear(); temp_deltas_concat.shrink_to_fit();    
        temp_bitmap_concat.clear(); temp_bitmap_concat.shrink_to_fit();
        temp_is_bitmap_marks.clear(); temp_is_bitmap_marks.shrink_to_fit();
        
    }

    int64_t serialize(ostream& os) const{
        int64_t bytes_written = 0;

        bytes_written += bitmap_concat.serialize(os);;
        bytes_written += bitmap_sizes.serialize(os);

        bytes_written += deltas_concat.serialize(os);;
        bytes_written += deltas_sizes.serialize(os);

        bytes_written += is_bitmap_marks.serialize(os);;
        bytes_written += is_bitmap_marks_rs.serialize(os);;

        return bytes_written;

        // Do not serialize temp structures
    }

    void load(istream& is){
        bitmap_concat.load(is);
        bitmap_sizes.load(is);
        deltas_concat.load(is);
        deltas_sizes.load(is);
        is_bitmap_marks.load(is);

        is_bitmap_marks_rs.load(is, &is_bitmap_marks);

        // Do not load temp structures
    }

    int64_t number_of_sets_stored() const{
        return is_bitmap_marks.size();
    }

    vector<New_Hybrid_Color_Set> get_all_sets() const{
        vector<New_Hybrid_Color_Set> all;
        for(int64_t i = 0; i < number_of_sets_stored(); i++){
            all.push_back(get_color_set_by_id(i));
        }
        return all;
    }


};