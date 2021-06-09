#pragma once

#include "stdint.h"
#include <cassert>
#include <bitset>
#include <unordered_map>

using namespace std;

// A bit-packed k-mer class for k <= 32.
class Kmer{
friend class std::hash<Kmer>;
friend class kmer_colex_compare;

private:
    uint64_t data; // Least significant bit pair is leftmost nucleotide
    uint8_t k;

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

    Kmer() : data(0), k(0) {}
    Kmer(uint8_t k) : data(0), k(k) {
        assert(k <= 32);
    }

    Kmer(uint64_t data, uint8_t k) : data(data), k(k) {
        assert(k <= 32);
    }

    Kmer(const string& S) : data(0), k(S.size()) {
        assert(S.size() <= 32);
        for(LL i = 0; i < (LL)S.size(); i++){
            set(i, S[i]);
        }
    }

    Kmer(const char* S, uint8_t k) : data(0), k(k) {
        assert(k <= 32);
        for(LL i = 0; i < k; i++){
            set(i, S[i]);
        }
    }

    uint8_t get_k() const{return k;}

    // Get the character at index idx from the start
    char get(LL idx) const{
        assert(idx >= 0 && idx < 32);
        if(idx >= k) return '\0';
        return to_char((data >> (idx*2)) & 0x03);
    }

    // Get the character at index idx from the start
    void set(LL idx, char c){
        assert(idx >= 0 && idx < 32);
        assert(c == 'A' || c == 'C' || c == 'G' || c == 'T');
        data &= (~(uint64_t)0) ^ ((uint64_t)0xff << (idx*2)); // Clear
        data |= (uint64_t)to_bitpair(c) << (idx*2); // Set
    }

    bool operator==(const Kmer &other) const {
        return this->data == other.data && this->k == other.k;
    }

    bool operator!=(const Kmer &other) const {
        return !(*this == other);
    }

    // Strict colexicographic comparison
    bool operator<(const Kmer& other) const{
        LL d = (int8_t)(this->k) - other.k;
        if(d == 0) return this->data < other.data; // Same k
        else if(d < 0){ // this is shorter
            if(this->data == (other.data >> abs(d)*2))
                return true; // this is a proper suffix of other
            else 
                return this->data < (other.data >> abs(d)*2);
        } else{ // other is shorter
            if((this->data >> abs(d)*2) == other.data)
                return false; // other is a proper suffix of this
            else 
                return (this->data >> abs(d)*2) < other.data;
        }
    }

    // Drops leftmost nucleotide, i.e. returns a (k-1)-mer (itself)
    Kmer dropleft(){
        assert(k > 0);
        data >>= 2;
        data &= (~(uint64_t)0) ^ ((uint64_t)0xff << ((k-1)*2)); // Clear
        k--;
        return *this;
    }

    // Drops rightmost nucleotide, i.e. returns a (k-1)-mer (itself)
    Kmer dropright(){
        assert(k > 0);
        data &= (~(uint64_t)0) ^ ((uint64_t)0xff << ((k-1)*2)); // Clear
        k--;
        return *this;
    }

    // Returns a k-mer like this but c added to the right (itself). Fails if
    // the length of the k-mer is 32.
    Kmer appendright(char c){
        assert(k < 32);
        k++;
        set(k-1, c);
        return *this;
    }

    // Returns a k-mer like this but c added to the left (itself). Fails if
    // the length of the k-mer is 32.
    Kmer appendleft(char c){
        assert(k < 32);
        k++;
        data <<= 2;
        set(0, c);
        return *this;
    }

    Kmer copy(){
        Kmer other(data, k);
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
        out.write((char*)(&data), sizeof(data));
        out.write((char*)(&k), sizeof(k));
    }

    void load(istream& in){
        in.read((char*)(&data), sizeof(data));
        in.read((char*)(&k), sizeof(k));
        assert(k <= 32);
    }

    void load(const char* bytes){
        memcpy(&data, bytes, sizeof(data));
        memcpy(&k, bytes + sizeof(data), sizeof(k));
        assert(k <= 32);
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
template<>
struct std::hash<Kmer>{
    std::size_t operator()(const Kmer& X) const {
        size_t h = std::hash<uint64_t>()(X.data);
        hash_combine<uint8_t>(h, X.k);
        return h;
    }
};

// Strict colexicographic comparison
struct kmer_colex_compare{
    inline bool operator() (const Kmer& A, const Kmer& B){
        return A < B;
    }
};