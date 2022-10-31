#pragma once

#include <vector>
#include "sdsl/bit_vectors.hpp"
#include "sdsl/int_vector.hpp"
#include "Color_Set_Interface.hh"
#include "SeqIO.hh"
#include <variant>

/*

This file defines a hybrid color set that is either a bit map or an integer array.

To avoid copying stuff around in memory, we have a color set view class that stores
just the pointer to the data. The pointer is stored as a std::variant, which stores 
either pointer and contains the bit specifying which type of pointer it is (bitmap or array).

But here is the problem. The user might want a mutable color set, for example when doing
intersections of the sets. The pointer will go into a static concatenation of sets
(see Color_Set_Storage.hh), which does not support creating new sets dynamically. This 
is why also we have a color set class that allocates and frees its own memory. This has 
exactly the same interface as the view class, except it has functions for set intersection 
and union, which modify the color set.

Now, to avoid duplicating code between the member functions of the color set view and the
concrete color set, we write the member functions as template functions that take the color
set (view or concrete) as the first parameter. The actual member functions in the color set
and view classes are just proxies that call these template functions. This saves 50+ lines
of duplicated code for member functions that are identical, such as get_colors_as_vector. 

You may be thinking, why not just use inheritance? That is not an option because the view
class needs a const-pointer whereas the other class needs a regular pointer. The view needs
a const-pointer so that the get-colorset method of the coloring storage can be const. The mutable
class needs a regular pointer because, well, it is supposed to be mutable. A std::variant of
a const-pointer and a regular pointer would be possible, but very ugly. Virtual inheritance
is not an option for performance reasons.

*/

using namespace std;

template<typename colorset_t> 
static inline bool colorset_is_empty(const colorset_t& cs){
    return cs.length == 0;
}

template<typename colorset_t> 
static inline bool colorset_is_bitmap(const colorset_t& cs){
    // Check variant type by index because it could be const or non-const and we don't want to care
    return cs.data_ptr.index() == 0;
}

template<typename colorset_t> 
static inline bool colorset_access_bitmap(const colorset_t& cs, int64_t idx){
    // Using std::holds_alternative by index because it could have a const or a non-const type
    return (*std::get<0>(cs.data_ptr))[cs.start + idx];
}

template<typename colorset_t> 
static inline int64_t colorset_access_array(const colorset_t& cs, int64_t idx){
    // Using std::holds_alternative by index because it could have a const or a non-const type
    return (*std::get<1>(cs.data_ptr))[cs.start + idx];
}

template<typename colorset_t> 
static inline int64_t colorset_size(const colorset_t& cs){
    if(colorset_is_bitmap(cs)){
        // Count number of bits set
        int64_t count = 0; 
        for(int64_t i = 0; i < cs.length; i++){
            count += colorset_access_bitmap(cs, i);
        }
        return count;
    } else return cs.length; // Array
}

template<typename colorset_t> 
static inline int64_t colorset_size_in_bits(const colorset_t& cs){
    if(colorset_is_bitmap(cs)) return cs.length;
    else return cs.length * std::get<1>(cs.data_ptr)->width();
    // Using std::holds_alternative by index because it could have a const or a non-const type
}

template<typename colorset_t> 
static inline vector<int64_t> colorset_get_colors_as_vector(const colorset_t& cs){
    std::vector<int64_t> vec;
    if(colorset_is_bitmap(cs)){    
        for(int64_t i = 0; i < cs.length; i++){
            if(colorset_access_bitmap(cs,i)) vec.push_back(i);
        }
    } else{
        for(int64_t i = 0; i < cs.length; i++){
            vec.push_back(colorset_access_array(cs,i));
        }
    }
    return vec;
}

template<typename colorset_t> 
static inline bool colorset_contains(const colorset_t& cs, int64_t color){
    if(colorset_is_bitmap(cs)){
        if(color >= cs.length) return false;
        return colorset_access_bitmap(cs, color);
    } else{
        // Linear scan. VERY SLOW
        for(int64_t i = 0; i < cs.length; i++){
            if(colorset_access_array(cs,i) == color) return true;
        }
        return false;
    }
}

// Stores the intersection into buf1 and returns the number of elements in the
// intersection (does not resize buf1). Buffer elements must be sorted.
// Assumes all elements in a buffer are distinct
int64_t intersect_buffers(sdsl::int_vector<>& buf1, int64_t buf1_len, const sdsl::int_vector<>& buf2, int64_t buf2_len);

// Stores the union into result_buf and returns the number of elements in the
// union (does not resize result_buf). Buffers elements must be sorted.
// Assumes all elements in a buffer are distinct. result_buf must have enough
// space to accommodate the union
int64_t union_buffers(vector<int64_t>& buf1, int64_t buf1_len, vector<int64_t>& buf2, int64_t buf2_len, vector<int64_t>& result_buf);

// Stores the result into A and returns the length of the new bit vector. A is not resized
// but the old elements past the end are left in place to avoid memory reallocations.
int64_t bitmap_vs_bitmap_intersection(sdsl::bit_vector& A, int64_t A_size, const sdsl::bit_vector& B, int64_t B_start, int64_t B_size);

// Stores the result into iv and returns the size of the intersection. iv is not resized
// but the old elements past the end  are left in place to avoid memory reallocations.
int64_t array_vs_bitmap_intersection(sdsl::int_vector<>& iv, int64_t iv_size, const sdsl::bit_vector& bv, int64_t bv_start, int64_t bv_size);

// Stores the result into bv and returns the length of bv. bv is not resized
// but the old elements past the end are left in place to avoid memory reallocations.
int64_t bitmap_vs_array_intersection(sdsl::bit_vector& bv, int64_t bv_size, const sdsl::int_vector<>& iv, int64_t iv_start, int64_t iv_size);

// Stores the result into A and returns the length of the new vector. A is not resized
// but the old elements past the end are left in place to avoid memory reallocations.
int64_t array_vs_array_intersection(sdsl::int_vector<>& A, int64_t A_len, const sdsl::int_vector<>& B, int64_t B_start, int64_t B_len);

// Stores the result into A and returns the length of the new bit vector. A must have enough
// space to accommodate the union
int64_t bitmap_vs_bitmap_union(sdsl::bit_vector& A, int64_t A_size, const sdsl::bit_vector& B, int64_t B_start, int64_t B_size);

// Stores the result into A and returns the length of the new bit vector. A must have enough
// space to accommodate the union
int64_t array_vs_bitmap_union(sdsl::int_vector<>& iv, int64_t iv_size, const sdsl::bit_vector& bv, int64_t bv_start, int64_t bv_size);

// Stores the result into A and returns the length of the new bit vector. A must have enough
// space to accommodate the union
int64_t bitmap_vs_array_union(sdsl::bit_vector& bv, int64_t bv_size, const sdsl::int_vector<>& iv, int64_t iv_start, int64_t iv_size);

// Stores the result into A and returns the length of the new bit vector. A must have enough
// space to accommodate the union
int64_t array_vs_array_union(sdsl::int_vector<>& A, int64_t A_len, const sdsl::int_vector<>& B, int64_t B_start, int64_t B_len);

class SDSL_Variant_Color_Set;

class SDSL_Variant_Color_Set_View{

public:

    std::variant<const sdsl::bit_vector*, const sdsl::int_vector<>*> data_ptr; // Non-owning pointer to external data
    int64_t start;
    int64_t length; // Number of bits in case of bit vector, number of elements in case of array

    SDSL_Variant_Color_Set_View(std::variant<const sdsl::bit_vector*, const sdsl::int_vector<>*> data_ptr, int64_t start, int64_t length)
        : data_ptr(data_ptr), start(start), length(length){}

    SDSL_Variant_Color_Set_View(const SDSL_Variant_Color_Set& cs); // Defined in the .cpp file because Color_Set is not yet defined at this point of this header

    bool empty() const {return colorset_is_empty(*this);};
    bool is_bitmap() const {return colorset_is_bitmap(*this);};
    int64_t size() const {return colorset_size(*this);}
    int64_t size_in_bits() const {return colorset_size_in_bits(*this);}
    bool contains(int64_t color) const {return colorset_contains(*this, color);}
    vector<int64_t> get_colors_as_vector() const {return colorset_get_colors_as_vector(*this);}

};

class SDSL_Variant_Color_Set{

    public:

    typedef SDSL_Variant_Color_Set_View view_t;

    std::variant<sdsl::bit_vector*, sdsl::int_vector<>*> data_ptr = (sdsl::bit_vector*) nullptr; // Owning pointer
    int64_t start = 0;
    int64_t length = 0; // Number of bits in case of bit vector, number of elements in case of array

    SDSL_Variant_Color_Set(){
        data_ptr = new sdsl::bit_vector();
    }

    const SDSL_Variant_Color_Set& operator=(const SDSL_Variant_Color_Set& other){
        if(this == &other) return *this; // Assignment to itself

        // Free existing memory
        auto call_delete = [](auto ptr){delete ptr;};
        std::visit(call_delete, this->data_ptr);

        // Alocate new memory
        if(std::holds_alternative<sdsl::bit_vector*>(other.data_ptr))
            this->data_ptr = new sdsl::bit_vector(*(std::get<sdsl::bit_vector*>(other.data_ptr)));
        else
            this->data_ptr = new sdsl::int_vector<>(*(std::get<sdsl::int_vector<>*>(other.data_ptr)));
        this->start = other.start;
        this->length = other.length;
        return *this;
    }

    SDSL_Variant_Color_Set(const SDSL_Variant_Color_Set& other){
        *this = other;
    }

    // Construct a copy of a color set from a view
    SDSL_Variant_Color_Set(const SDSL_Variant_Color_Set_View& view) : length(view.length){
        if(std::holds_alternative<const sdsl::bit_vector*>(view.data_ptr)){
            data_ptr = new sdsl::bit_vector(view.length, 0);

            // Copy the bits
            const sdsl::bit_vector* from = std::get<const sdsl::bit_vector*>(view.data_ptr);
            sdsl::bit_vector* to = std::get<sdsl::bit_vector*>(data_ptr);
            for(int64_t i = 0; i < view.length; i++){
                (*to)[i] = (*from)[view.start + i];
            }
        } else{
            // Array
            int64_t bit_width = std::get<const sdsl::int_vector<>*>(view.data_ptr)->width();
            data_ptr = new sdsl::int_vector<>(view.length, 0, bit_width);

            // Copy the values
            const sdsl::int_vector<>* from = std::get<const sdsl::int_vector<>*>(view.data_ptr);
            sdsl::int_vector<>* to = std::get<sdsl::int_vector<>*>(data_ptr);
            for(int64_t i = 0; i < view.length; i++){
                (*to)[i] = (*from)[view.start + i];
            }
        }
    }

    SDSL_Variant_Color_Set(const vector<int64_t>& set){
        int64_t max_element = *std::max_element(set.begin(), set.end());
        if(log2(max_element) * set.size() > max_element){ // TODO: this if-statement is duplicated in this file
            // Dense -> bitmap
            sdsl::bit_vector* ptr = new sdsl::bit_vector(max_element+1, 0);
            for(int64_t x : set) (*ptr)[x] = 1;
            length = ptr->size();
            data_ptr = ptr; // Assign to variant
        } else{
            // Sparse -> array
            sdsl::int_vector<>* ptr = new sdsl::int_vector<>(set.size());
            if(set.size() > 0){
                for(int64_t i = 0; i < set.size(); i++){
                    (*ptr)[i] = set[i];
                }
            }
            length = ptr->size();
            data_ptr = ptr; // Assign to variant
        }
    }

    ~SDSL_Variant_Color_Set(){
        auto call_delete = [](auto ptr){delete ptr;};
        std::visit(call_delete, data_ptr);
    }

    bool empty() const {return colorset_is_empty(*this);};
    bool is_bitmap() const {return colorset_is_bitmap(*this);};
    int64_t size() const {return colorset_size(*this);}
    int64_t size_in_bits() const {return colorset_size_in_bits(*this);}
    bool contains(int64_t color) const {return colorset_contains(*this, color);}
    vector<int64_t> get_colors_as_vector() const {return colorset_get_colors_as_vector(*this);}

    // Stores the intersection back to to this object
    void intersection(const SDSL_Variant_Color_Set_View& other){
        if(is_bitmap() && other.is_bitmap()){
            this->length = bitmap_vs_bitmap_intersection(*std::get<sdsl::bit_vector*>(data_ptr), this->length, *std::get<const sdsl::bit_vector*>(other.data_ptr), other.start, other.length);
        } else if(!is_bitmap() && other.is_bitmap()){
            this->length = array_vs_bitmap_intersection(*std::get<sdsl::int_vector<>*>(data_ptr), this->length, *std::get<const sdsl::bit_vector*>(other.data_ptr), other.start, other.length);
        } else if(is_bitmap() && !other.is_bitmap()){
            // The result will be sparse, so this will turn our representation into an array
            sdsl::int_vector<> iv_copy(*std::get<const sdsl::int_vector<>*>(other.data_ptr)); // Make a mutable copy

            // Intersect into the mutable copy
            int64_t iv_copy_length = array_vs_bitmap_intersection(iv_copy, other.length, *std::get<sdsl::bit_vector*>(this->data_ptr), this->start, this->length);

            // Replace our data with the mutable copy
            auto call_delete = [](auto ptr){delete ptr;}; // Free current data
            std::visit(call_delete, this->data_ptr);
            this->data_ptr = new sdsl::int_vector<>(iv_copy);
            this->length = iv_copy_length;
        } else{ // Array vs Array
            this->length = array_vs_array_intersection(*std::get<sdsl::int_vector<>*>(data_ptr), this->length, *std::get<const sdsl::int_vector<>*>(other.data_ptr), other.start, other.length);
        }
    }

    void do_union(const SDSL_Variant_Color_Set_View& other){
        // TODO: DO PROPERLY
        vector<int64_t> A = this->get_colors_as_vector();
        vector<int64_t> B = other.get_colors_as_vector();
        vector<int64_t> AB(A.size() + B.size());
        int64_t len = union_buffers(A, A.size(), B, B.size(), AB);
        AB.resize(len);
        *this = SDSL_Variant_Color_Set(AB);
    }

};
