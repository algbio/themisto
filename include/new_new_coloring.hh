#pragma once

#include <vector>
#include "sdsl/bit_vectors.hpp"
#include "sdsl/int_vector.hpp"
#include <variant>

using namespace std;

// This class either owns its memory, or points to other memory
class New_Hybrid_Color_Set{
    
private:

    bool access_bitmap(int64_t idx) const{
        return (*std::get<sdsl::bit_vector*>(data_ptr))[start + idx];
    }

    int64_t access_delta_array(int64_t idx) const{
        return (*std::get<sdsl::int_vector<>*>(data_ptr))[start + idx];
    }

public:

    int64_t start;
    int64_t length; // Number of bits in case of bit vector, number of elements in case of delta array

    bool owns_memory;
    std::variant<sdsl::bit_vector*, sdsl::int_vector<>*> data_ptr; // Pointer to external data

    bool is_bitmap() const{
        return std::holds_alternative<sdsl::bit_vector*>(data_ptr);
    }

    // Construct with pointer to outside memory
    New_Hybrid_Color_Set(int64_t start, int64_t length, std::variant<sdsl::bit_vector*, sdsl::int_vector<>*> data_ptr)
        : start(start), length(length), owns_memory(false), data_ptr(data_ptr) {}

    // Construct with own memory. Always a bitmap for now TODO
    // v must be sorted
    New_Hybrid_Color_Set(const vector<int64_t>& v) : start(0), length(v.back() + 1), owns_memory(true){
        data_ptr = new sdsl::bit_vector(length, 0);
        for(int64_t x : v) (*std::get<sdsl::bit_vector*>(data_ptr))[x] = 1;
    }

    // TODO: need copy construction and assignment operator because of memory management. Otherwise might
    // end up with double delete.

    New_Hybrid_Color_Set() : start(0), length(0), owns_memory(false){}

    ~New_Hybrid_Color_Set(){
        if(owns_memory){
            if(std::holds_alternative<sdsl::bit_vector*>(data_ptr))
                delete std::get<sdsl::bit_vector*>(data_ptr);
            else
                delete std::get<sdsl::int_vector<>*>(data_ptr);
        }
    }

    bool empty() const{
        return length == 0;
    }

    int64_t size() const{
        if(is_bitmap()){
            // Count number of bits set
            int64_t count = 0; 
            for(int64_t i = 0; i < length; i++){
                count += access_bitmap(i);
            }
            return count;
        } else return length; // Array
    }

    int64_t size_in_bits() const{
        if(is_bitmap()) return length;
        else return length * std::get<sdsl::int_vector<>*>(data_ptr)->width();
    }

    bool contains(int64_t color) const{
        if(is_bitmap()){
            if(color >= length) return false;
            return access_bitmap(color);
        } else{
            // Linear scan. VERY SLOWS
            for(int64_t x : get_colors_as_vector()) if(x == color) return true;
            return false;
        }
    }

    New_Hybrid_Color_Set intersection(){
        return New_Hybrid_Color_Set(); // TODO
    }

    New_Hybrid_Color_Set do_union(){
        return New_Hybrid_Color_Set(); // TODO
    }

    int64_t serialize() const{
        return 0; // TODO
    }

    void load() const{
        return; //TODO
    }

    vector<int64_t> get_colors_as_vector() const{
        std::vector<int64_t> vec;
        if(is_bitmap()){    
            for(int64_t i = 0; i < length; i++){
                if(access_bitmap(i)) vec.push_back(i);
            }
        } else{
            for(int64_t i = 0; i < length; i++){
                if(i == 0) vec.push_back(access_delta_array(0));
                else vec.push_back(vec.back() + access_delta_array(i));
            }
        }
        return vec;
    }

};

template<typename color_set_t> class Color_Set_Storage;

template<>
class Color_Set_Storage<New_Hybrid_Color_Set>{

    private:

    sdsl::bit_vector bitmap_concat;
    sdsl::int_vector<> bitmap_starts; // bitmap_starts[i] = starting position of the i-th bitmap

    sdsl::int_vector<> deltas_concat;
    sdsl::int_vector<> deltas_starts; // deltas_starts[i] = starting position of the i-th delta list

    sdsl::bit_vector is_bitmap_marks;
    sdsl::rank_support_v5<> is_bitmap_marks_rs;

    // Dynamic-length vectors used during construction only
    // TODO: refactor these out of the class to a separate construction class
    vector<bool> temp_bitmap_concat;
    vector<int64_t> temp_deltas_concat;
    vector<int64_t> temp_bitmap_starts;
    vector<int64_t> temp_deltas_starts;
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

    sdsl::int_vector<> to_sdsl_int_vector(const vector<int64_t>& v){
        int64_t max_element = *std::max_element(v.begin(), v.end());
        sdsl::int_vector iv(v.size(), 0, bits_needed(max_element));
        for(int64_t i = 0; i < v.size(); i++) iv[i] = v[i];
        return iv;
    }

    public:

    New_Hybrid_Color_Set get_color_set_by_id(int64_t id) const{
        if(is_bitmap_marks[id]){
            int64_t bitmap_idx = is_bitmap_marks_rs.rank(id); // This many bitmaps come before this bitmap
            int64_t start = bitmap_starts[bitmap_idx];
            int64_t end = bitmap_starts[bitmap_idx+1]; // One past the end

            // The following const cast is dangerous but it is ok because the color set class will not modify
            // the pointer. The cast is needed because the color set class also manages its possible own memory
            // at the same pointer, which it needs to free.
            std::variant<sdsl::bit_vector*, sdsl::int_vector<>*> data_ptr = const_cast<sdsl::bit_vector*>(&bitmap_concat);
            return New_Hybrid_Color_Set(start, end-start, data_ptr);
        } else{
            int64_t deltas_idx = id - is_bitmap_marks_rs.rank(id); // Rank-0. This many delta arrays come before this bitmap
            int64_t start = deltas_starts[deltas_idx];
            int64_t end = deltas_starts[deltas_idx+1]; // One past the end

            // See comment about const-cast above
            std::variant<sdsl::bit_vector*, sdsl::int_vector<>*> data_ptr = const_cast<sdsl::int_vector<>*>(&deltas_concat);
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

            // Store bitmap start
            temp_bitmap_starts.push_back(temp_bitmap_concat.size());

            // Create bitmap
            vector<bool> bitmap(set.size(), 0);
            for(int64_t x : set) bitmap[x] = 1;
            for(bool b : bitmap) temp_bitmap_concat.push_back(b);

        } else{
            // Sparse -> delta array

            // Add is_bitmap_mark
            temp_is_bitmap_marks.push_back(0);

            // Store array start
            temp_deltas_starts.push_back(temp_deltas_concat.size());

            if(set.size() > 0){
                for(int64_t i = 0; i < set.size(); i++){
                    if(i == 0) temp_deltas_concat.push_back(set[0]);
                    else temp_deltas_concat.push_back(set[i] - set[i-1]);
                }
                
            }
        }

    }


    // Call this after done with add_set
    void prepare_for_queries(){

        // Add extra starts points one past the end
        // These eliminate a special case when querying for the size of the last color set
        temp_bitmap_starts.push_back(temp_bitmap_starts.size());
        temp_deltas_starts.push_back(temp_deltas_starts.size());

        deltas_concat = to_sdsl_int_vector(temp_deltas_concat);
        bitmap_starts = to_sdsl_int_vector(temp_bitmap_starts);
        deltas_starts = to_sdsl_int_vector(temp_deltas_starts);
        bitmap_concat = to_sdsl_bit_vector(temp_bitmap_concat);
        is_bitmap_marks = to_sdsl_bit_vector(temp_is_bitmap_marks);

        sdsl::util::init_support(is_bitmap_marks_rs, &is_bitmap_marks);

        // Free memory
        temp_deltas_concat.clear(); temp_deltas_concat.shrink_to_fit();    
        temp_bitmap_concat.clear(); temp_bitmap_concat.shrink_to_fit();
        temp_is_bitmap_marks.clear(); temp_is_bitmap_marks.shrink_to_fit();
        temp_deltas_starts.clear(); temp_deltas_starts.shrink_to_fit();
        temp_bitmap_starts.clear(); temp_bitmap_starts.shrink_to_fit();
    }

    int64_t serialize(ostream& os) const{
        int64_t bytes_written = 0;

        bytes_written += bitmap_concat.serialize(os);
        bytes_written += bitmap_starts.serialize(os);

        bytes_written += deltas_concat.serialize(os);
        bytes_written += deltas_starts.serialize(os);

        bytes_written += is_bitmap_marks.serialize(os);
        bytes_written += is_bitmap_marks_rs.serialize(os);

        return bytes_written;

        // Do not serialize temp structures
    }

    void load(istream& is){
        bitmap_concat.load(is);
        bitmap_starts.load(is);
        deltas_concat.load(is);
        deltas_starts.load(is);
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