#pragma once

#include "stdint.h"
#include <string.h> // memcpy
#include <cassert>
#include <bitset>
#include <unordered_map>
#include <iostream>
#include <vector>

using namespace std;

// Table mapping ACGT to 0123
static constexpr uint8_t from_ACGT_to_0123_lookup_table[256] = 
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// Table mapping ACGT to 0123
static constexpr char from_0123_to_ACGT_lookup_table[4] = {'A','C','G','T'};

template<int64_t max_len> class kmer_colex_compare; // Compiler needs this forward declaration to make it a friend class

// A bit-packed k-mer class containing at most max_len characters
// Maximum length supported is 255 because the length is stored in an uint8_t
template<int64_t max_len> 
class Kmer{
friend class std::hash<Kmer<max_len>>;
friend class kmer_colex_compare<max_len>;

private:

    // Two bits per character.
    // The rightmost character of the whole k-mer is at the most significant bits of block 0.
    // See the function `get` for a precise definition of the data layout.
    static constexpr int64_t DATA_ARRAY_SIZE = max_len /  32 + (max_len % 32 > 0);
    uint64_t data[DATA_ARRAY_SIZE]; 
    uint8_t k;
    
    char constexpr to_char(uint8_t x) const{
        assert(x < 4);
        return from_0123_to_ACGT_lookup_table[x]; // In globals.hh
    }

    char constexpr to_bitpair(char c) const{
        return from_ACGT_to_0123_lookup_table[(uint8_t)c];  // In globals.hh
    }

public:

    Kmer() : k(0) {
        clear();
    }
    
    Kmer(uint8_t k) : k(k) {
        assert(k <= min(max_len, (int64_t)255));
        clear();
    }

    Kmer(const string& S) : k(S.size()) {
        assert(S.size() <= min(max_len, (int64_t)255));
        clear();
        for(int64_t i = 0; i < (int64_t)S.size(); i++){
            set(i, S[i]);
        }
    }

    Kmer(const char* S, uint8_t k) : k(k) {
        assert(k <= min(max_len, (int64_t)255));
        clear();
        for(int64_t i = 0; i < k; i++){
            set(i, S[i]);
        }
    }

    void clear(){
        for(int64_t i = 0; i < DATA_ARRAY_SIZE; i++) data[i] = 0;
    }

    uint8_t get_k() const{return k;}

    // Get the character at index idx from the start
    char get(int64_t idx) const{
        assert(idx < max_len);
        if(idx >= k) return '\0';
        int64_t block_idx = (k-1-idx) / 32;
        int64_t block_offset = (k-1-idx) % 32;
        return to_char((data[block_idx] >> ((31 - block_offset)*2)) & 0x03);
    }

    // Get the character at index idx from the start
    void set(int64_t idx, char c){
        assert(idx >= 0 && idx < max_len);
        assert(c == 'A' || c == 'C' || c == 'G' || c == 'T');
        int64_t block_idx = (k-1-idx) / 32;
        int64_t block_offset = (k-1-idx) % 32;
        data[block_idx] &= (~(uint64_t)0) ^ ((uint64_t)0x03 << ((31 - block_offset)*2)); // Clear
        data[block_idx] |= (uint64_t)to_bitpair(c) << ((31 - block_offset)*2); // Set
    }

    bool operator==(const Kmer &other) const {
        if(this->k != other.k) return false;
        for(int64_t i = 0; i < DATA_ARRAY_SIZE; i++)
            if((this->data)[i] != other.data[i]) return false;
        return true;
    }

    bool operator!=(const Kmer &other) const {
        return !(*this == other);
    }

    // Strict colexicographic comparison
    // Returns true if this is smaller than the other
    bool operator<(const Kmer& other) const{
        // If there is a mismatch in the data blocks, then
        // case 1: the mismatch is in the k-mer part of both k-mers. Then the mismatch determines the order.
        // case 2: the mismatch is past the end of the shorter k-mer. The shorter k-mer has an implicit 'A'
        // in that position, so the longer k-mer must has something larger than 'A'. This decides the mismatch
        // in favor of the shorter k-mer, which is correct.
        for(int64_t i = 0; i < DATA_ARRAY_SIZE; i++){
            if(this->data[i] < other.data[i]) return true;
            if(this->data[i] > other.data[i]) return false;
        }
        // All data blocks were equal. If the k-mer lengths are equal, then the k-mers are equal, so the answer is false.
        // If the k-mer lengths are not equal, then one is a prefix of the othe (and the part of the other that is longer
        // is all A-characters). Then, the shorter k-mer is the smaller.
        if(this->k < other.k) return true;
        else return false;
    }

    bool operator>(const Kmer& other) const{
        return !(*this < other) && !(*this == other);
    }

    // Drops leftmost nucleotide, i.e. returns a (k-1)-mer (itself)
    Kmer dropleft(){
        assert(k > 0);
        // Shift block
        int64_t block_idx = (k-1) / 32;
        int64_t block_offset = (k-1) % 32;
        data[block_idx] &= (~(uint64_t)0) ^ ((uint64_t)0x03 << ((31 - block_offset)*2)); // Clear
        k--;
        return *this;
    }

    // Drops rightmost nucleotide, i.e. returns a (k-1)-mer (itself)
    Kmer dropright(){
        assert(k > 0);
        for(int64_t block = 0; block < DATA_ARRAY_SIZE; block++){
            data[block] <<= 2; // Shift
            // Carry character from the next block
            if(block != DATA_ARRAY_SIZE-1)
                data[block] |= (data[block+1] >> (31*2)) & (uint64_t)0x03;
        }
        k--;
        return *this;
    }

    // Returns a k-mer like this but c added to the right (itself). Fails if
    // the length of the k-mer is 32.
    Kmer appendright(char c){
        assert(k < max_len);
        k++;
        // Shift blocks to make room
        for(int64_t i = DATA_ARRAY_SIZE-1; i >= 0; i--){
            data[i] >>= 2;
            if(i != 0){
                // Carry last character of the previous block
                data[i] |= (data[i-1] & (uint64_t)0x03) << (31*2);
            }
        }
        set(k-1, c);
        return *this;
    }

    // Returns a k-mer like this but c added to the left (itself). Fails if
    // the length of the k-mer is 32.
    Kmer appendleft(char c){
        assert(k < max_len);
        k++;
        set(0, c); // append c to the left
        return *this;
    }

    Kmer copy() const{
        Kmer other;
        other.k = this->k;
        for(int64_t i = 0; i < DATA_ARRAY_SIZE; i++)
            other.data[i] = (this->data)[i];
        return other;
    }

    string to_string() const{
        string S;
        for(int64_t i = 0; i < k; i++) S += get(i);
        return S;
    }

    char last() const{
        return get(get_k()-1);
    }

    char first() const{
        return get(0);
    }

    static int64_t size_in_bytes(){
        return sizeof(data) + sizeof(k);
    }

    void serialize(ostream& out) const{
        out.write((char*)(data), sizeof(data));
        out.write((char*)(&k), sizeof(k));
    }

    void load(istream& in){
        in.read((char*)(data), sizeof(data));
        in.read((char*)(&k), sizeof(k));
        assert(k <= min(max_len, (int64_t)255));
    }


    // `out` must have at least size_in_bytes() bytes of space
    void serialize(char* out) const{
        memcpy(out, (char*)(data), sizeof(data));
        memcpy(out + sizeof(data), (char*)(&k), sizeof(k));
    }

    // Load from memory serialized to by member function `serialize(char* out)`
    void load(const char* in){
        memcpy((char*)(data), in, sizeof(data));
        memcpy((char*)(&k), in + sizeof(data), sizeof(k));
        assert(k <= min(max_len, (int64_t)255));
    }

};

// https://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x
template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

// Strict colexicographic comparison
template<int64_t max_len>
struct kmer_colex_compare{
    inline bool operator() (const Kmer<max_len>& A, const Kmer<max_len>& B){
        return A < B;
    }
};

// For unordered_map
template<int64_t max_len>
struct std::hash<Kmer<max_len>>{
    std::size_t operator()(const Kmer<max_len>& X) const {
        size_t h = 0;
        for(int64_t i = 0; i < Kmer<max_len>::DATA_ARRAY_SIZE; i++)
            hash_combine<uint64_t>(h, X.data[i]);
        hash_combine<uint8_t>(h, X.k);
        return h;
    }
};
