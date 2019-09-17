#ifndef SDSL_ITERATE_HH
#define SDSL_ITERATE_HH

#include <sdsl/construct.hpp>
#include <vector>
#include <utility>
#include <string>
#include "bwt.hh"
#include "Interval.hh"
#include "sdsl/io.hpp"
#include "sdsl/wavelet_trees.hpp"
//#include "../../Interfaces.hh"


/*
 * Implements a bidictional BWT index for a byte alphabet.
 * All indices and ranks are indexed starting from zero.
 * All intervals are inclusive, i.e. the interval from i to j includes both i and j.
 *
 * Terminology: C-array of an interval is an array with size 256 such that
 * C[c] is the number of characters with lexicographical rank strictly less than c at the given interval. 
 * If the interval spans the whole BWT, the array is called the global C-array,
 * else it's called a local C-array.
 */

template<class t_bitvector = sdsl::bit_vector>
class BD_BWT_index{

public:

    struct Interval_Data {

        sdsl::int_vector_size_type n_distinct_symbols;
        std::vector<uint8_t> symbols; // Vector of length n_distinct_symbols

        // ranks_start[i] = number of symbols smaller than the i-th symbol before the start of the interval
        std::vector<uint64_t> ranks_start;

        // ranks_start[i] = number of symbols smaller than the i-th symbol up to the end of the interval
        std::vector<uint64_t> ranks_end;

        int64_t count(int64_t symbol_index){
            return ranks_end[symbol_index] - ranks_start[symbol_index];
        }

    };


private:
        
    sdsl::wt_hutu<t_bitvector> forward_bwt;
    sdsl::wt_hutu<t_bitvector> reverse_bwt;

    std::vector<int64_t> global_c_array;
    std::vector<uint8_t> alphabet;
    
    // Space reused between iterations
    std::vector<int64_t> local_c_array;
    std::vector<uint8_t> symbols;
    std::vector<uint64_t> ranks_i;
    std::vector<uint64_t> ranks_j;
    
    std::vector<uint8_t> get_string_alphabet(const uint8_t* s) const;
    int64_t strlen(const uint8_t* str) const;
    int64_t compute_cumulative_char_rank_in_interval(const sdsl::wt_hutu<t_bitvector>& wt, uint8_t c, Interval I) const;
    int64_t get_interval_symbols(const sdsl::wt_hutu<t_bitvector>& wt, Interval I);
    void get_interval_symbols(const sdsl::wt_hutu<t_bitvector>& wt, Interval I, sdsl::int_vector_size_type& nExtensions, 
                              std::vector<uint8_t>& symbols, std::vector<uint64_t>& ranks_i, std::vector<uint64_t>& ranks_j) const;
    void count_smaller_chars(const sdsl::wt_hutu<t_bitvector>& bwt, std::vector<int64_t>& counts, Interval I);

public:

    // The input string must not contain the END byte
    static const uint8_t END = 0x01; // End of string marker.
    BD_BWT_index() {}
    BD_BWT_index(const uint8_t* input);
    
    uint8_t get_END() const { return END; }
    int64_t size() const { return forward_bwt.size();}
    uint8_t forward_bwt_at(int64_t index) const { return forward_bwt[index]; }
    uint8_t backward_bwt_at(int64_t index) const { return reverse_bwt[index]; }
    const std::vector<int64_t>& get_global_c_array() const { return global_c_array; }
    const std::vector<uint8_t>& get_alphabet() const { return alphabet; }

    // Computes the local C-array of the given forward interval into the parameter vector. The size
    // of the parameter vector must be at least 256
    void compute_local_c_array_forward(Interval& interval, std::vector<int64_t>& c_array);

    // Computes the local C-array of the given reverse interval into the parameter vector. The size
    // of the parameter vector must be at least 256
    void compute_local_c_array_reverse(Interval& interval, std::vector<int64_t>& c_array);

    // Computes the interval pair of the left extension of the given intervals by the character c.
    // Returns an interval of size zero if the extension is not possible or if the given interval has size 0.
    // Internally, the function computes the local C-array of the interval. 
    // To avoid recomputing the local C-array multiple times for the same interval, 
    // use the version of the function that takes the C-array as a parameter.
    Interval_pair left_extend(Interval_pair intervals, uint8_t c);
    
    // A version of left_extend that takes a precomputed local forward c-array as a parameter. Useful
    // to avoid recomputing the same local c-array for the same interval.
    Interval_pair left_extend(Interval_pair intervals, uint8_t c, const std::vector<int64_t>& local_c_array);

    // Analogous right extensions to left extensions
    Interval_pair right_extend(Interval_pair intervals, uint8_t c);
    // Takes a precomputed local reverse c-array as a parameter.
    Interval_pair right_extend(Interval_pair intervals, uint8_t c, const std::vector<int64_t>& local_c_array);

    // Let s be a suffix of length k with lexicographic rank lex_rank among all suffixes of the input string.
    // Returns the lexicographic rank of the suffix with length k+1 if k is not equal to the length of the
    // input string counting the terminating symbol, or the lexicographic rank of the suffix with length 0 otherwise.
    int64_t backward_step(int64_t lex_rank) const;
    
    // Let s be a prefix of length k with colexicographic rank colex_rank among all prefixes of the input string.
    // Returns the colexicographical rank of the prefix with length k+1 if k is not equal to the length of the
    // input string counting the terminating symbol, or the colexicographic rank of the prefix with length 0 otherwise.
    int64_t forward_step(int64_t colex_rank) const;
    
    bool is_right_maximal(Interval_pair I);
    bool is_left_maximal(Interval_pair I);

    // Results are stored in the provided struct reference
    void compute_bwt_interval_data(Interval I, Interval_Data& data){
        if(I.size() == 0)
            data.n_distinct_symbols = 0;
        else 
            forward_bwt.interval_symbols(I.left, I.right+1, data.n_distinct_symbols, data.symbols, data.ranks_start, data.ranks_end);
            //       wt.interval_symbols(I.left, I.right+1, nExtensions, symbols, ranks_i, ranks_j);
    }
    
    // Results are stored in the provided struct reference
    void compute_rev_bwt_interval_data(Interval I, Interval_Data& data){
        if(I.size() == 0)
            data.n_distinct_symbols = 0;
        else 
            reverse_bwt.interval_symbols(I.left, I.right+1, data.n_distinct_symbols, data.symbols, data.ranks_start, data.ranks_end);
    }
    
    // Symbol index: Position of the extensions symbols in data.symbols
    Interval_pair left_extend(Interval_pair intervals, Interval_Data& data, int64_t symbol_index){
        if(intervals.forward.size() == 0)
            return Interval_pair(-1,-2,-1,-2);
        
        char c = data.symbols[symbol_index];
        
        // Compute the new forward interval
        int64_t num_c_in_interval = data.ranks_end[symbol_index] - data.ranks_start[symbol_index];
        
        int64_t start_f_new = get_global_c_array()[c] + data.ranks_start[symbol_index]; // Start in forward
        int64_t end_f_new = start_f_new + num_c_in_interval - 1; // End in forward

        if(start_f_new > end_f_new) return Interval_pair(-1,-2,-1,-2); // num_c_in_interval == 0
        
        // Compute the new reverse interval
        Interval reverse = intervals.reverse;
        int64_t smaller_chars = 0;
        for(int64_t i = 0; i < symbol_index; i++){
            smaller_chars += data.ranks_end[i] - data.ranks_start[i];
        }
        
        int64_t start_r_new = reverse.left + smaller_chars;
        int64_t end_r_new = start_r_new + (end_f_new - start_f_new); // The forward and reverse intervals must have same length
        
        return Interval_pair(start_f_new,end_f_new,start_r_new,end_r_new);
    }
    
    Interval_pair right_extend(Interval_pair intervals, Interval_Data& data, int64_t symbol_index){
        if(intervals.reverse.size() == 0)
            return Interval_pair(-1,-2,-1,-2);
        
        char c = data.symbols[symbol_index];
        
        // Compute the new reverse interval
        int64_t num_c_in_interval = data.ranks_end[symbol_index] - data.ranks_start[symbol_index];
        
        int64_t start_r_new = get_global_c_array()[c] + data.ranks_start[symbol_index]; // Start in forward
        int64_t end_r_new = start_r_new + num_c_in_interval - 1; // End in forward

        if(start_r_new > end_r_new) return Interval_pair(-1,-2,-1,-2); // num_c_in_interval == 0
        
        // Compute the new forward interval
        Interval forward = intervals.forward;
        int64_t smaller_chars = 0;
        for(int64_t i = 0; i < symbol_index; i++){
            smaller_chars += data.ranks_end[i] - data.ranks_start[i];
        }
        
        int64_t start_f_new = forward.left + smaller_chars;
        int64_t end_f_new = start_f_new + (end_r_new - start_r_new); // The forward and reverse intervals must have same length
        
        return Interval_pair(start_f_new,end_f_new,start_r_new,end_r_new);
    }
    
    void save_to_disk(std::string directory, std::string filename_prefix);
    void save_to_disk_reverse_only(std::string directory, std::string filename_prefix);
    void load_from_disk(std::string directory, std::string filename_prefix);

};

template<class t_bitvector>
const uint8_t BD_BWT_index<t_bitvector>::END;

template<class t_bitvector>
void BD_BWT_index<t_bitvector>::compute_local_c_array_forward(Interval& interval, std::vector<int64_t>& c_array){
    assert(c_array.size() >= 256);
    count_smaller_chars(forward_bwt, c_array, interval);
}

template<class t_bitvector>
void BD_BWT_index<t_bitvector>::compute_local_c_array_reverse(Interval& interval, std::vector<int64_t>& c_array){
    assert(c_array.size() >= 256);
    count_smaller_chars(reverse_bwt, c_array, interval);
}

template<class t_bitvector>
int64_t BD_BWT_index<t_bitvector>::compute_cumulative_char_rank_in_interval(const sdsl::wt_hutu<t_bitvector>& wt, uint8_t c, Interval I) const{
    int64_t ans = 0;
    if(I.size() == 0) return 0;
    
    // Sum of ranks of all characters that are lexicographically smaller than c
    for(uint8_t d : alphabet){
        if(d == c) break;
        ans += wt.rank(I.right + 1,d) - wt.rank(I.left,d);
    }
    return ans;
}


// Returns the number of distinct symbols and stores them in this->symbols without resizising this->symbols
template<class t_bitvector>
int64_t BD_BWT_index<t_bitvector>::get_interval_symbols(const sdsl::wt_hutu<t_bitvector>& wt, Interval I){
    if(I.size() == 0){
        std::vector<uint8_t> empty;
        return 0;
    }
    
    sdsl::int_vector_size_type nExtensions;
    wt.interval_symbols(I.left, I.right+1, nExtensions, this->symbols, this->ranks_i, this->ranks_j);
    return nExtensions;
}

// Another version of get_interval_symbols when the allocation of a new vector is too slow.
// Counts the number of distinct symbols into nExtensions, and puts the symbols into the indices
// [0,nExtensions[ of symbols. Also stores ranks of the symbols at the endpoints of the interval I
// to ranks_i and ranks_j. Important: All the parameter vectors must have length at least equal to the size
// of the alphabet of the given wavelet tree. Symbols is not sorted
template<class t_bitvector>
void BD_BWT_index<t_bitvector>::get_interval_symbols(const sdsl::wt_hutu<t_bitvector>& wt, Interval I, sdsl::int_vector_size_type& nExtensions, std::vector<uint8_t>& symbols,
 std::vector<uint64_t>& ranks_i, std::vector<uint64_t>& ranks_j) const{
    if(I.size() == 0){
        nExtensions = 0;
        return;
    }
    
    wt.interval_symbols(I.left, I.right+1, nExtensions, symbols, ranks_i, ranks_j);
}


// Takes a backward step in the forward bwt
template<class t_bitvector>
int64_t BD_BWT_index<t_bitvector>::backward_step(int64_t lex_rank) const{
    uint8_t c = forward_bwt[lex_rank];
    return global_c_array[c] + forward_bwt.rank(lex_rank, c);
}

// Takes a backward step in the reverse bwt
template<class t_bitvector>
int64_t BD_BWT_index<t_bitvector>::forward_step(int64_t colex_rank) const{
    uint8_t c = reverse_bwt[colex_rank];
    return global_c_array[c] + reverse_bwt.rank(colex_rank, c);
}

template<class t_bitvector>
Interval_pair BD_BWT_index<t_bitvector>::left_extend(Interval_pair intervals, uint8_t c){
    compute_local_c_array_forward(intervals.forward, local_c_array);
    return left_extend(intervals,c,local_c_array);
}


template<class t_bitvector>
Interval_pair BD_BWT_index<t_bitvector>::right_extend(Interval_pair intervals, uint8_t c){
    compute_local_c_array_reverse(intervals.reverse, local_c_array);
    return right_extend(intervals,c,local_c_array);
}

template<class t_bitvector>
Interval_pair BD_BWT_index<t_bitvector>::left_extend(Interval_pair intervals, uint8_t c, const std::vector<int64_t>& local_c_array){
    assert(local_c_array.size() >= 256);
    if(intervals.forward.size() == 0)
        return Interval_pair(-1,-2,-1,-2);
    
    Interval forward = intervals.forward;
    Interval reverse = intervals.reverse;
    
    // Compute the new forward interval
    int64_t num_c_in_interval = forward_bwt.rank(forward.right + 1,c) - forward_bwt.rank(forward.left,c);
    int64_t start_f_new = get_global_c_array()[c] + forward_bwt.rank(forward.left, c); // Start in forward
    int64_t end_f_new = start_f_new + num_c_in_interval - 1; // End in forward

    if(start_f_new > end_f_new) return Interval_pair(-1,-2,-1,-2); // num_c_in_interval == 0
    
    // Compute the new reverse interval
    int64_t start_r_new = reverse.left + local_c_array[c];
    int64_t end_r_new = start_r_new + (end_f_new - start_f_new); // The forward and reverse intervals must have same length
    
    return Interval_pair(start_f_new,end_f_new,start_r_new,end_r_new);
}

template<class t_bitvector>
Interval_pair BD_BWT_index<t_bitvector>::right_extend(Interval_pair intervals, uint8_t c, const std::vector<int64_t>& local_c_array){
    assert(local_c_array.size() >= 256);
    if(intervals.forward.size() == 0)
        return Interval_pair(-1,-2,-1,-2);
    
    Interval forward = intervals.forward;
    Interval reverse = intervals.reverse;
    
    // Compute the new reverse interval
    int64_t num_c_in_interval = reverse_bwt.rank(reverse.right + 1,c) - reverse_bwt.rank(reverse.left,c);
    int64_t start_r_new = get_global_c_array()[c] + reverse_bwt.rank(reverse.left, c); // Start in reverse
    int64_t end_r_new = start_r_new + num_c_in_interval - 1; // End in reverse

    if(start_r_new > end_r_new) return Interval_pair(-1,-2,-1,-2); // num_c_in_interval == 0
    
    // Compute the new forward interval
    int64_t start_f_new = forward.left + local_c_array[c];
    int64_t end_f_new = start_f_new + (end_r_new - start_r_new); // The forward and reverse intervals must have same length
    
    return Interval_pair(start_f_new,end_f_new,start_r_new,end_r_new);
    
}

// counts = vector with 256 elements
// Modifies the counts vector such that:
// counts[i] = number of characters in interval I that are lexicographically stricly smaller than alphabet[i]
// Assumes alphabet is sorted

template<class t_bitvector>
void BD_BWT_index<t_bitvector>::count_smaller_chars(const sdsl::wt_hutu<t_bitvector>& bwt, 
                                                    std::vector<int64_t>& counts, Interval I){
    assert(alphabet.size() != 0);
    for(int64_t i = 0; i < (int64_t)alphabet.size(); i++) counts[alphabet[i]] = 0;
    if(I.size() == 0) return;
    
    sdsl::int_vector_size_type nSymbols;
    bwt.interval_symbols(I.left, I.right+1, nSymbols, symbols, ranks_i, ranks_j);
    
    // Put in the counts
    for(int64_t k = 0; k < nSymbols; k++){
        counts[symbols[k]] = ranks_j[k] - ranks_i[k];
    }    
    
    // Do the cumulative sum
    int64_t prev_cumul = 0;
    int64_t prev_count = counts[alphabet[0]];
    int64_t temp = 0;
    for(int64_t k = 0; k < alphabet.size(); k++){
        temp = counts[alphabet[k]];
        if(k == 0) counts[alphabet[k]] = 0;
        else{
            counts[alphabet[k]] = prev_cumul + prev_count;
        }
        prev_count = temp;
        prev_cumul = counts[alphabet[k]];
    }    
}

template<class t_bitvector>
bool BD_BWT_index<t_bitvector>::is_right_maximal(Interval_pair I){
    
    // An interval is right-maximal iff it has more than one possible right extension
    int64_t nSymbols = get_interval_symbols(reverse_bwt, I.reverse);
    return nSymbols >= 2;
}

template<class t_bitvector>
bool BD_BWT_index<t_bitvector>::is_left_maximal(Interval_pair I){
    
    // An interval is left-maximal iff it has more than one possible left extension
    int64_t nSymbols = get_interval_symbols(forward_bwt, I.forward);
    return nSymbols >= 2;
}

// Returns the alphabet in sorted order
template<class t_bitvector>
std::vector<uint8_t> BD_BWT_index<t_bitvector>::get_string_alphabet(const uint8_t* s) const{
    
    std::vector<bool> found(256,false);
    while(*s != 0){
        found[*s] = true;
        s++;
    }

    std::vector<uint8_t> alphabet;
    for(int i = 0; i < 256; i++){
        if(found[i]) alphabet.push_back((uint8_t)i);
    }
    
    return alphabet;
}

// strlen(const uint8_t*) is not in the standard library
template<class t_bitvector>
int64_t BD_BWT_index<t_bitvector>::strlen(const uint8_t* str) const{
    const uint8_t* start = str;
    while(*str != 0) str++;
    return str - start;
}

template<class t_bitvector>
BD_BWT_index<t_bitvector>::BD_BWT_index(const uint8_t* input) 
    : global_c_array(256), local_c_array(256), symbols(256), ranks_i(256), ranks_j(256) {
    if(*input == 0) throw std::runtime_error("Tried to construct BD_BWT_index for an empty string");
    int64_t n = strlen(input);
    
    if(std::find(input, input+n, END) != input + n){
        std::stringstream error;
        error << "Input string contains forbidden byte " << std::hex << END;
        throw std::runtime_error(error.str());
    }
    
    // Build the two bwts
    uint8_t* forward = (uint8_t*) malloc(sizeof(uint8_t) * (n + 1));
    uint8_t* backward = (uint8_t*) malloc(sizeof(uint8_t) * (n + 1));
    for(int64_t i = 0; i < n; i++){
        forward[i] = input[i];
        backward[n-1-i] = input[i];
    }
    forward[n] = END;
    backward[n] = END;

    uint8_t* forward_transform = build_bwt(forward,n,END);
    free(forward);
        
    uint8_t* backward_transform = build_bwt(backward,n,END);
    free(backward);
    
    // Build wavelet trees
    construct_im(this->forward_bwt, (const char*)forward_transform, 1); // Must cast to signed char* or else breaks. File a bug report to sdsl?
    construct_im(this->reverse_bwt, (const char*)backward_transform, 1); // Must cast to signed char* or else breaks. File a bug report to sdsl?
    
    this->alphabet = get_string_alphabet(forward_transform);
    
    free(forward_transform);
    free(backward_transform); 
    
    // Compute cumulative character counts
    count_smaller_chars(forward_bwt,global_c_array,Interval(0,forward_bwt.size()-1));
}

template<class t_bitvector>
void BD_BWT_index<t_bitvector>::save_to_disk(std::string directory, std::string filename_prefix){
    std::string fbwt = directory + "/" + filename_prefix + "_forward_bwt.dat";
    if(!sdsl::store_to_file(forward_bwt, fbwt)) {
        throw std::runtime_error("Error writing to disk: " + fbwt);
    }
    
    std::string rbwt = directory + "/" + filename_prefix + "_reverse_bwt.dat";
    if(!sdsl::store_to_file(reverse_bwt, rbwt)) {
        throw std::runtime_error("Error writing to disk: " + rbwt);
    }
    
    // Copy to sdsl bit vector because they have serialization built in
    sdsl::int_vector<64> global_c_array_sdsl(global_c_array.size());
    for(int64_t i = 0; i < global_c_array.size(); i++){
        global_c_array_sdsl[i] = global_c_array[i];
    }
    std::string gca = directory + "/" + filename_prefix + "_gca.dat";
    if(!sdsl::store_to_file(global_c_array_sdsl, gca)) {
        throw std::runtime_error("Error writing to disk: " + gca);
    }
        
    // Copy to sdsl bit vector because they have serialization built in
    sdsl::int_vector<8> alphabet_sdsl(alphabet.size());
    for(int64_t i = 0; i < alphabet.size(); i++){
        alphabet_sdsl[i] = alphabet[i];
    }
    std::string A = directory + "/" + filename_prefix + "_alphabet.dat";
    if(!sdsl::store_to_file(alphabet_sdsl, A)) {
        throw std::runtime_error("Error writing to disk: " + A);
    }
}

template<class t_bitvector>
void BD_BWT_index<t_bitvector>::save_to_disk_reverse_only(std::string directory, std::string filename_prefix){
    
    std::string rbwt = directory + "/" + filename_prefix + "_reverse_bwt.dat";
    if(!sdsl::store_to_file(reverse_bwt, rbwt)) {
        throw std::runtime_error("Error writing to disk: " + rbwt);
    }
    
    // Copy to sdsl bit vector because they have serialization built in
    sdsl::int_vector<64> global_c_array_sdsl(global_c_array.size());
    for(int64_t i = 0; i < global_c_array.size(); i++){
        global_c_array_sdsl[i] = global_c_array[i];
    }
    std::string gca = directory + "/" + filename_prefix + "_gca.dat";
    if(!sdsl::store_to_file(global_c_array_sdsl, gca)) {
        throw std::runtime_error("Error writing to disk: " + gca);
    }
        
    // Copy to sdsl bit vector because they have serialization built in
    sdsl::int_vector<8> alphabet_sdsl(alphabet.size());
    for(int64_t i = 0; i < alphabet.size(); i++){
        alphabet_sdsl[i] = alphabet[i];
    }
    std::string A = directory + "/" + filename_prefix + "_alphabet.dat";
    if(!sdsl::store_to_file(alphabet_sdsl, A)) {
        throw std::runtime_error("Error writing to disk: " + A);
    }
}


template<class t_bitvector>
void BD_BWT_index<t_bitvector>::load_from_disk(std::string directory, std::string filename_prefix){
    global_c_array.resize(256); 
    local_c_array.resize(256); 
    symbols.resize(256); 
    ranks_i.resize(256);
    ranks_j.resize(256);
    
    std::string fbwt = directory + "/" + filename_prefix + "_forward_bwt.dat";
    if(!sdsl::load_from_file(forward_bwt, fbwt)) {
        throw std::runtime_error("Error reading from disk: " + fbwt);
    }
    std::string rbwt = directory + "/" + filename_prefix + "_reverse_bwt.dat";
    if(!sdsl::load_from_file(reverse_bwt, rbwt)) {
        throw std::runtime_error("Error reading from disk: " + rbwt);
    }
    
    sdsl::int_vector<64> gca_sdsl;
    std::string gca_path = directory + "/" + filename_prefix + "_gca.dat";
    if(!sdsl::load_from_file(gca_sdsl, gca_path)) {
        throw std::runtime_error("Error reading from disk: " + gca_path);
    }
    for(int64_t i = 0; i < 256; i++){
        global_c_array[i] = gca_sdsl[i];
    }
    
    sdsl::int_vector<8> alphabet_sdsl;
    std::string alphabet_path = directory + "/" + filename_prefix + "_alphabet.dat";
    if(!sdsl::load_from_file(alphabet_sdsl, alphabet_path)) {
        throw std::runtime_error("Error reading from disk: " + alphabet_path);
    }
    
    alphabet.resize(alphabet_sdsl.size());
    for(int64_t i = 0; i < alphabet_sdsl.size(); i++){
        alphabet[i] = alphabet_sdsl[i];
    }
}

#endif
