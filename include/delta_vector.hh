#pragma once

#include <vector>
#include <sdsl/coder_elias_delta.hpp>

using namespace std;

class Fixed_Width_Delta_Vector{

    // Number of bits required to represent x
    int64_t bits_needed(uint64_t x){
        return max((int64_t)std::bit_width(x), (int64_t)1); // Need at least 1 bit (for zero)
    }


public:

    sdsl::int_vector<> diffs; // Diffs data[i] = v[i] - v[i-1], such that data[0] := v[0]

    Delta_Vector() {}

    // The input values must be sorted and distinct
    Delta_Vector(const vector<int64_t>& sorted_values){
        if(sorted_values.size() > 0){
            int64_t largest_diff = sorted_values[0];
            for(int64_t i = 1; i < sorted_values.size(); i++)
                largest_diff = std::max(largest_diff, sorted_values[i] - sorted_values[i-1]);

            diffs = sdsl::int_vector<>(sorted_values.size(), 0, bits_needed(largest_diff));

            diffs[0] = sorted_values[0];
            for(int64_t i = 1; i < sorted_values.size(); i++)
                diffs[i] = sorted_values[i] - sorted_values[i-1];
        }
    }

    vector<int64_t> get_values() const{
        vector<int64_t> values(diffs.size());
        for(int64_t i = 0; i < values.size(); i++){
            if(i == 0) values[i] = diffs[0];
            else values[i] = values[i-1] + diffs[i];
        }
        return values;
    }

    bool empty() const {
        return diffs.empty();
    }

    int64_t size_in_bytes() const{
        return sdsl::size_in_bytes(diffs);
    }

    int64_t serialize(std::ostream& os) const{
        int64_t n_bytes_written = 0;

        n_bytes_written += diffs.serialize(os);

        return n_bytes_written;
    }

    void load(std::istream& is){
        diffs.load(is);
    }

};

class Elias_Delta_Vector{

public:

    sdsl::int_vector<> data; // Elias-delta-encoded diffs data[i] = v[i] - v[i-1], such that data[0] := v[0] + 1

    Delta_Vector() {}

    // The input values must be sorted and distinct
    Delta_Vector(const vector<int64_t>& sorted_values){
        if(sorted_values.size() > 0){
            sdsl::int_vector<> diffs(sorted_values.size());
            diffs[0] = sorted_values[0] + 1; // Elias-delta does not have a code for zero, so we add one
            for(int64_t i = 1; i < sorted_values.size(); i++)
                diffs[i] = sorted_values[i] - sorted_values[i-1];
            sdsl::coder::elias_delta::encode(diffs, data);
        }
    }

    vector<int64_t> get_values() const{
        sdsl::int_vector<> diffs;
        sdsl::coder::elias_delta::decode(data, diffs);
        vector<int64_t> values(diffs.size());
        for(int64_t i = 0; i < values.size(); i++){
            if(i == 0) values[i] = diffs[0] - 1; // Undo the added +1 in the constructor
            else values[i] = values[i-1] + diffs[i];
        }
        return values;
    }

    bool empty() const {
        return data.bit_size() == 0;
    }

    int64_t size_in_bytes() const{
        return sdsl::size_in_bytes(data);
    }

    int64_t serialize(std::ostream& os) const{
        int64_t n_bytes_written = 0;

        n_bytes_written += data.serialize(os);

        return n_bytes_written;
    }

    void load(std::istream& is){
        data.load(is);
    }

};