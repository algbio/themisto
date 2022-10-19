#pragma once

#include <vector>
#include "sdsl/bit_vectors.hpp"
#include "sdsl/int_vector.hpp"

using namespace std;

class New_Hybrid_Color_Set{
    bool is_bitmap;
    
    int64_t bitmap_start;
    int64_t bitmap_length;

    int64_t deltas_start;
    int64_t deltas_end;

    bool owns_memory;

    sdsl::bit_vector* bv_ptr;
    sdsl::int_vector<>* iv_ptr;

};

template<typename color_set_t> class Color_Set_Storage;

template<>
class Color_Set_Storage<New_Hybrid_Color_Set>{

    private:

    sdsl::bit_vector bitmap_concat;
    sdsl::bit_vector bitmap_unary_sizes; // In units of bits
    sdsl::rank_support_v5<> bitmap_unary_sizes_rs;
    sdsl::select_support_mcl<> bitmap_unary_sizes_ss;

    sdsl::int_vector<> deltas_concat;
    sdsl::bit_vector deltas_unary_sizes; // In units of elements
    sdsl::rank_support_v5<> deltas_unary_sizes_rs;
    sdsl::select_support_mcl<> deltas_unary_sizes_ss;

    sdsl::bit_vector is_bitmap_marks;
    sdsl::rank_support_v5<> is_bitmap_marks_rs;

    // Dynamic-length vectors used during construction only
    vector<bool> temp_bitmap_concat;
    vector<bool> temp_bitmap_unary_sizes;

    vector<int64_t> temp_deltas_concat;
    vector<bool> temp_deltas_unary_sizes;

    vector<bool> temp_is_bitmap_marks;

    void append_unary_number(int64_t x, vector<bool>& v){
        for(int64_t i = 0; i < x; i++) v.push_back(0);
        v.push_back(1);
    }

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

    const New_Hybrid_Color_Set& get_color_set_by_id(int64_t id) const{
        return {};
        // TODO
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
            append_unary_number(set.size(), temp_bitmap_unary_sizes);
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
            append_unary_number(set.size(), temp_deltas_unary_sizes);
        }

    }


    // Call this after done with add_set
    void prepare_for_queries(){
        int64_t max_delta = *std::max_element(temp_deltas_concat.begin(), temp_deltas_concat.end());
        deltas_concat = sdsl::int_vector<>(temp_deltas_concat.size(), 0, bits_needed(max_delta));
        for(int64_t i = 0; i < temp_deltas_concat.size(); i++)
            deltas_concat[i] = temp_deltas_concat[i];

        bitmap_concat = to_sdsl_bit_vector(temp_bitmap_concat);
        bitmap_unary_sizes = to_sdsl_bit_vector(temp_bitmap_unary_sizes);
        deltas_unary_sizes = to_sdsl_bit_vector(temp_deltas_unary_sizes);
        is_bitmap_marks = to_sdsl_bit_vector(temp_is_bitmap_marks);

        sdsl::util::init_support(bitmap_unary_sizes_rs, &bitmap_unary_sizes);
        sdsl::util::init_support(bitmap_unary_sizes_ss, &bitmap_unary_sizes);
        sdsl::util::init_support(deltas_unary_sizes_rs, &deltas_unary_sizes);
        sdsl::util::init_support(deltas_unary_sizes_ss, &deltas_unary_sizes);
        sdsl::util::init_support(is_bitmap_marks_rs, &is_bitmap_marks);

        // Free memory
        temp_deltas_concat.clear(); temp_deltas_concat.shrink_to_fit();    
        temp_bitmap_concat.clear(); temp_bitmap_concat.shrink_to_fit();
        temp_bitmap_unary_sizes.clear(); temp_bitmap_unary_sizes.shrink_to_fit();
        temp_deltas_unary_sizes.clear(); temp_deltas_unary_sizes.shrink_to_fit();
        temp_is_bitmap_marks.clear(); temp_is_bitmap_marks.shrink_to_fit();
        
    }

    int64_t serialize(ostream& os) const{
        int64_t byte_written = 0;

        bytes_written += bitmap_concat.serialize(os);;
        bytes_written += bitmap_unary_sizes.serialize(os);; // In units of bits
        bytes_written += bitmap_unary_sizes_rs.serialize(os);;
        bytes_written += bitmap_unary_sizes_ss.serialize(os);;

        bytes_written += deltas_concat.serialize(os);;
        bytes_written += deltas_unary_sizes.serialize(os);; // In units of elements
        bytes_written += deltas_unary_sizes_rs.serialize(os);;
        bytes_written += deltas_unary_sizes_ss.serialize(os);;

        bytes_written += is_bitmap_marks.serialize(os);;
        bytes_written += is_bitmap_marks_rs.serialize(os);;

        return bytes_written;

        // Do not serialize temp structures
    }

    void load(istream& is){
        bitmap_concat.load(is);
        bitmap_unary_sizes.load(is);
        bitmap_unary_sizes_rs.load(is);
        bitmap_unary_sizes_ss.load(is);
        deltas_concat.load(is);
        deltas_unary_sizes.load(is);
        deltas_unary_sizes_rs.load(is);
        deltas_unary_sizes_ss.load(is);
        is_bitmap_marks.load(is);
        is_bitmap_marks_rs.load(is);

        // Do not load temp structures
    }

    int64_t number_of_sets_stored() const{
        return is_bitmap_marks.size();
    }

    vector<New_Hybrid_Color_Set> get_all_sets() const{
        vector<New_Hybrid_Color_Set> all;
        for(int64_t i = 0; i < number_of_sets_stored; i++){
            all.push_back(get_color_set_by_id(i));
        }
        return all;
    }


};