#pragma once

#include <vector>
#include "sdsl/bit_vectors.hpp"
#include "sdsl/int_vector.hpp"

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
