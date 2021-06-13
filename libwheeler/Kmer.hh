#pragma once

#include "stdint.h"
#include "../globals.hh"
#include <string.h> // memcpy
#include <cassert>
#include <bitset>
#include <unordered_map>
#include <iostream>
#include <vector>

using namespace std;

template<int64_t max_len> class kmer_colex_compare; // Compiler needs this forward declaration to make it a friend class

// A bit-packed k-mer class containing at most max_len characters
// MAXIMUM LENGTH SUPPORTED IS 255 because the length is stored in an uint8_t
template<int64_t max_len> 
class Kmer{
friend class std::hash<Kmer<max_len>>;
friend class kmer_colex_compare<max_len>;

private:

    // Two bits per character. Least significant bit pair is leftmost nucleotide.
    // Like this: If 0 is the leftmost index in the string and n-1 is the rightmost,
    // then the characters are arranged in the the blocks like this (for illustration block have size 3 here):
    //
    // block 0 | block 1 | block 2
    //  8 7 6     5 4 3     2 1 0
    //
    // With this arrangement, colex comparison of equal-length k-mers
    // is the lex-comparison of the blocks from left to right, such that
    // each block is considered a single character. This means we can bit-
    // parallelize the colex comparison.
    static const int64_t DATA_ARRAY_SIZE = max_len /  32 + (max_len % 32 > 0);
    uint64_t data[DATA_ARRAY_SIZE]; 
    uint8_t k; // this is an 8-bit value to save space so the maximum k is 2^8 - 1 = 255

    char to_char(uint8_t x) const{
        switch(x){
            case 0x00: return 'A';
            case 0x01: return 'C';
            case 0x02: return 'G';
            case 0x03: return 'T';
        }
        assert(false);
        return 0;
    }

    char to_bitpair(char c) const{
        switch(c){
            case 'A': return 0x00;
            case 'C': return 0x01;
            case 'G': return 0x02;
            case 'T': return 0x03;
        }
        assert(false);
        return 0;
    }

public:

    Kmer() : k(0) {
        clear();
    }
    Kmer(uint8_t k) : k(k) {
        assert(k <= min(max_len, (LL)255));
        clear();
    }

    Kmer(const string& S) : k(S.size()) {
        assert(S.size() <= min(max_len, (LL)255));
        clear();
        for(LL i = 0; i < (LL)S.size(); i++){
            set(i, S[i]);
        }
    }

    Kmer(const char* S, uint8_t k) : k(k) {
        assert(k <= min(max_len, (LL)255));
        clear();
        for(LL i = 0; i < k; i++){
            set(i, S[i]);
        }
    }

    void clear(){
        for(LL i = 0; i < DATA_ARRAY_SIZE; i++) data[i] = 0;
    }

    uint8_t get_k() const{return k;}

    // Get the character at index idx from the start
    char get(LL idx) const{
        assert(idx < max_len);
        if(idx >= k) return '\0';
        LL block_idx = DATA_ARRAY_SIZE - 1 - idx / 32;
        LL block_offset = idx % 32;
        return to_char((data[block_idx] >> (block_offset*2)) & 0x03);
    }

    // Get the character at index idx from the start
    void set(LL idx, char c){
        assert(idx >= 0 && idx < max_len);
        assert(c == 'A' || c == 'C' || c == 'G' || c == 'T');
        LL block_idx = DATA_ARRAY_SIZE - 1 - idx / 32;
        LL block_offset = idx % 32;
        data[block_idx] &= (~(uint64_t)0) ^ ((uint64_t)0x03 << (block_offset*2)); // Clear
        data[block_idx] |= (uint64_t)to_bitpair(c) << (block_offset*2); // Set
    }

    bool operator==(const Kmer &other) const {
        if(this->k != other.k) return false;
        for(LL i = 0; i < DATA_ARRAY_SIZE; i++)
            if((this->data)[i] != other.data[i]) return false;
        return true;
    }

    bool operator!=(const Kmer &other) const {
        return !(*this == other);
    }

    // Strict colexicographic comparison
    bool operator<(const Kmer& other) const{
        LL n = min(this->k, other.k);
        for(LL i = 0; i < n; i++){
            char c1 = this->get(this->k - 1 - i);
            char c2 = other.get(other.k - 1 - i);
            if(c1 < c2) return true;
            else if(c1 > c2) return false;
            else continue;
        }
        // One of the k-mers is a prefix of the other. Return the shorter
        if(this->k < other.k) return true;
        else return false;
    }

    // Drops leftmost nucleotide, i.e. returns a (k-1)-mer (itself)
    Kmer dropleft(){
        assert(k > 0);
        // Shift block
        for(LL block = DATA_ARRAY_SIZE - 1; block >= 0; block--){
            data[block] >>= 2;
            if(block != 0){ // Carry a character from the next block
                data[block] |= (data[block-1] & 0x03) << (31*2); // Set
            }
        }

        // Clear the leftmost character
        LL last_block = DATA_ARRAY_SIZE - 1 - (k-1) / 32;
        LL last_char_offset = (k-1) % 32;
        data[last_block] &= (~(uint64_t)0) ^ ((uint64_t)0x03 << (last_char_offset*2)); // Clear
        k--;
        return *this;
    }

    // Drops rightmost nucleotide, i.e. returns a (k-1)-mer (itself)
    Kmer dropright(){
        assert(k > 0);
        LL last_block = DATA_ARRAY_SIZE - 1 - (k-1) / 32;
        LL last_char_offset = (k-1) % 32;
        data[last_block] &= (~(uint64_t)0) ^ ((uint64_t)0x03 << (last_char_offset*2)); // Clear
        k--;
        return *this;
    }

    // Returns a k-mer like this but c added to the right (itself). Fails if
    // the length of the k-mer is 32.
    Kmer appendright(char c){
        assert(k < max_len);
        k++;
        set(k-1, c);
        return *this;
    }

    // Returns a k-mer like this but c added to the left (itself). Fails if
    // the length of the k-mer is 32.
    Kmer appendleft(char c){
        assert(k < max_len);
        k++;
        for(LL i = 0; i < DATA_ARRAY_SIZE; i++){
            data[i] <<= 2;
            if(i != DATA_ARRAY_SIZE-1){
                // Carry last character of the next block
                data[i] |= (data[i+1] >> (31*2)) & (uint64_t)0x03;
            }
        }
        set(0, c); // append c to the left
        return *this;
    }

    Kmer copy(){
        Kmer other;
        other.k = this->k;
        for(LL i = 0; i < DATA_ARRAY_SIZE; i++)
            other.data[i] = (this->data)[i];
        return other;
    }

    string to_string() const{
        string S;
        for(LL i = 0; i < k; i++) S += get(i);
        return S;
    }

    char last() const{
        return get(get_k()-1);
    }

    char first() const{
        return get(0);
    }

    static LL size_in_bytes(){
        return sizeof(data) + sizeof(k);
    }

    void serialize(ostream& out) const{
        out.write((char*)(data), sizeof(data));
        out.write((char*)(&k), sizeof(k));
    }

    void load(istream& in){
        in.read((char*)(data), sizeof(data));
        in.read((char*)(&k), sizeof(k));
        assert(k <= min(max_len, (LL)255));
    }

    void load(const char* bytes){
        memcpy(data, bytes, sizeof(data));
        memcpy(&k, bytes + sizeof(data), sizeof(k));
        assert(k <= min(max_len, (LL)255));
    }

};

// https://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x
template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

// For unordered_map
template<int64_t max_len>
struct std::hash<Kmer<max_len>>{
    std::size_t operator()(const Kmer<max_len>& X) const {
        size_t h = std::hash<uint64_t>()(X.data);
        hash_combine<uint8_t>(h, X.k);
        return h;
    }
};

// Strict colexicographic comparison
template<int64_t max_len>
struct kmer_colex_compare{
    inline bool operator() (const Kmer<max_len>& A, const Kmer<max_len>& B){
        return A < B;
    }
};