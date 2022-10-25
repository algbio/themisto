#pragma once

#include <vector>
#include "Color_Set_Interface.hh"
#include "Color_Set.hh"
#include <iostream>
#include <map>

using namespace std;

// Takes as parameter a class that encodes a single color set, and a viewer class for that
template<typename colorset_t, typename colorset_view_t> 
requires Color_Set_Interface<colorset_t>
class Color_Set_Storage{

private:

    vector<colorset_t> sets;

public:

    Color_Set_Storage(){}
    Color_Set_Storage(const vector<colorset_t>& sets) : sets(sets) {
        prepare_for_queries();
    }

    colorset_view_t get_color_set_by_id(int64_t id) const{
        return colorset_view_t(sets[id]);
    }

    // Need to call prepare_for_queries() after all sets have been added
    void add_set(const vector<int64_t>& set){
        sets.push_back(set);
    }

    // Call this after done with add_set
    void prepare_for_queries(){
        sets.shrink_to_fit();
    }

    int64_t serialize(ostream& os) const{
        int64_t bytes_written = 0;
        std::size_t n_sets = sets.size();
        os.write(reinterpret_cast<char*>(&n_sets), sizeof(std::size_t));
        bytes_written += sizeof(std::size_t);

        for (std::size_t i = 0; i < n_sets; ++i) {
            bytes_written += sets[i].serialize(os);
        }
        return bytes_written;
    }

    void load(istream& is){
        std::size_t n_sets = 0;
        is.read(reinterpret_cast<char*>(&n_sets), sizeof(std::size_t));

        sets.resize(n_sets);
        for (std::size_t i = 0; i < n_sets; ++i) {
            colorset_t cs;
            cs.load(is);
            sets[i] = cs;
        }
    }

    int64_t number_of_sets_stored() const{
        return sets.size();
    }

    vector<colorset_view_t> get_all_sets() const{
        return sets;
    }

    // Returns map: component -> number of bytes
    map<string, int64_t> space_breakdown() const{
        map<string, int64_t> breakdown;
        sbwt::SeqIO::NullStream ns;
        int64_t total_set_byte_size = 0;
        for(int64_t i = 0; i < sets.size(); i++){
            total_set_byte_size += sets[i].serialize(ns);
        }
        breakdown["sets"] = total_set_byte_size;

        return breakdown;
    }
};


/* 

Template specialization for Color_Set and Color_Set_View

The purpose of this storage class is to have a color set class that avoids space
overheads associated with allocating each color set from the heap separately.
Instead, here we concatenate all the color sets and store pointers to the starts
of the color sets.

Actually, there are two concatenations: a concatenation of bit maps and a concatenation 
of integer arrays. A color is stored either as a bit map or an integer array, depending
which one is smaller.

*/

template<>
class Color_Set_Storage<Color_Set, Color_Set_View>{

    private:

    sdsl::bit_vector bitmap_concat;
    sdsl::int_vector<> bitmap_starts; // bitmap_starts[i] = starting position of the i-th bitmap

    sdsl::int_vector<> arrays_concat;
    sdsl::int_vector<> arrays_starts; // arrays_starts[i] = starting position of the i-th subarray

    sdsl::bit_vector is_bitmap_marks;
    sdsl::rank_support_v5<> is_bitmap_marks_rs;

    // Dynamic-length vectors used during construction only
    // TODO: refactor these out of the class to a separate construction class
    vector<bool> temp_bitmap_concat;
    vector<int64_t> temp_arrays_concat;
    vector<int64_t> temp_bitmap_starts;
    vector<int64_t> temp_arrays_starts;
    vector<bool> temp_is_bitmap_marks;

    // Number of bits required to represent x
    int64_t bits_needed(uint64_t x){
        return max((int64_t)std::bit_width(x), (int64_t)1); // Need at least 1 bit (for zero)
    }

    sdsl::bit_vector to_sdsl_bit_vector(const vector<bool>& v){
        if(v.size() == 0) return sdsl::bit_vector();
        sdsl::bit_vector bv(v.size());
        for(int64_t i = 0; i < v.size(); i++) bv[i] = v[i];
        return bv;
    }

    sdsl::int_vector<> to_sdsl_int_vector(const vector<int64_t>& v){
        if(v.size() == 0) return sdsl::int_vector<>();
        int64_t max_element = *std::max_element(v.begin(), v.end());
        sdsl::int_vector iv(v.size(), 0, bits_needed(max_element));
        for(int64_t i = 0; i < v.size(); i++) iv[i] = v[i];
        return iv;
    }

    public:

    Color_Set_Storage() {}

    // Build from list of sets directly. Instead of using this constructor, you should probably just
    // call add_set for each set you want to add separately and then call prepare_for_queries when done.
    // This constructor is just to have the same interface as the other color set storage class.
    Color_Set_Storage(const vector<Color_Set>& sets){
        for(const Color_Set& cs : sets){
            add_set(cs.get_colors_as_vector());
        }
        prepare_for_queries();
    }

    Color_Set_View get_color_set_by_id(int64_t id) const{
        if(is_bitmap_marks[id]){
            int64_t bitmap_idx = is_bitmap_marks_rs.rank(id); // This many bitmaps come before this bitmap
            int64_t start = bitmap_starts[bitmap_idx];
            int64_t end = bitmap_starts[bitmap_idx+1]; // One past the end

            std::variant<const sdsl::bit_vector*, const sdsl::int_vector<>*> data_ptr = &bitmap_concat;
            return Color_Set_View(data_ptr, start, end-start);
        } else{
            int64_t arrays_idx = id - is_bitmap_marks_rs.rank(id); // Rank-0. This many arrays come before this bitmap
            int64_t start = arrays_starts[arrays_idx];
            int64_t end = arrays_starts[arrays_idx+1]; // One past the end

            std::variant<const sdsl::bit_vector*, const sdsl::int_vector<>*> data_ptr = &arrays_concat;
            return Color_Set_View(data_ptr, start, end-start);
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
            vector<bool> bitmap(max_element+1, 0);
            for(int64_t x : set) bitmap[x] = 1;
            for(bool b : bitmap) temp_bitmap_concat.push_back(b);

        } else{
            // Sparse -> Array

            // Add is_bitmap_mark
            temp_is_bitmap_marks.push_back(0);

            // Store array start
            temp_arrays_starts.push_back(temp_arrays_concat.size());

            if(set.size() > 0){
                for(int64_t i = 0; i < set.size(); i++){
                    temp_arrays_concat.push_back(set[i]);
                }
            }
        }

    }


    // Call this after done with add_set
    void prepare_for_queries(){

        // Add extra starts points one past the end
        // These eliminate a special case when querying for the size of the last color set
        temp_bitmap_starts.push_back(temp_bitmap_concat.size());
        temp_arrays_starts.push_back(temp_arrays_concat.size());

        arrays_concat = to_sdsl_int_vector(temp_arrays_concat);
        bitmap_starts = to_sdsl_int_vector(temp_bitmap_starts);
        arrays_starts = to_sdsl_int_vector(temp_arrays_starts);
        bitmap_concat = to_sdsl_bit_vector(temp_bitmap_concat);
        is_bitmap_marks = to_sdsl_bit_vector(temp_is_bitmap_marks);

        sdsl::util::init_support(is_bitmap_marks_rs, &is_bitmap_marks);

        // Free memory
        temp_arrays_concat.clear(); temp_arrays_concat.shrink_to_fit();    
        temp_bitmap_concat.clear(); temp_bitmap_concat.shrink_to_fit();
        temp_is_bitmap_marks.clear(); temp_is_bitmap_marks.shrink_to_fit();
        temp_arrays_starts.clear(); temp_arrays_starts.shrink_to_fit();
        temp_bitmap_starts.clear(); temp_bitmap_starts.shrink_to_fit();
    }

    int64_t serialize(ostream& os) const{
        int64_t bytes_written = 0;

        bytes_written += bitmap_concat.serialize(os);
        bytes_written += bitmap_starts.serialize(os);

        bytes_written += arrays_concat.serialize(os);
        bytes_written += arrays_starts.serialize(os);

        bytes_written += is_bitmap_marks.serialize(os);
        bytes_written += is_bitmap_marks_rs.serialize(os);

        return bytes_written;

        // Do not serialize temp structures
    }

    void load(istream& is){
        bitmap_concat.load(is);
        bitmap_starts.load(is);
        arrays_concat.load(is);
        arrays_starts.load(is);
        is_bitmap_marks.load(is);

        is_bitmap_marks_rs.load(is, &is_bitmap_marks);

        // Do not load temp structures
    }

    int64_t number_of_sets_stored() const{
        return is_bitmap_marks.size();
    }

    vector<Color_Set_View> get_all_sets() const{
        vector<Color_Set_View> all;
        for(int64_t i = 0; i < number_of_sets_stored(); i++){
            all.push_back(get_color_set_by_id(i));
        }
        return all;
    }

    // Returns map: component -> number of bytes
    map<string, int64_t> space_breakdown() const{
        map<string, int64_t> breakdown;

        sbwt::SeqIO::NullStream ns;

        breakdown["bitmaps-concat"] = bitmap_concat.serialize(ns);
        breakdown["bitmaps-starts"] = bitmap_starts.serialize(ns);
        breakdown["arrays-concat"] = arrays_concat.serialize(ns);
        breakdown["arrays-starts"] = arrays_starts.serialize(ns);
        breakdown["is-bitmap-marks"] = is_bitmap_marks.serialize(ns);
        breakdown["is-bitmap-marks-rank-suppport"] = is_bitmap_marks_rs.serialize(ns);

        // In the future maybe the space breakdown struct should support float statistics but for now we just disgustingly print to cout.
        cout << "Fraction of bitmaps in coloring: " << (double) is_bitmap_marks_rs.rank(is_bitmap_marks.size()) / is_bitmap_marks.size() << endl;

        return breakdown;
    }

};