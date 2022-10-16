#pragma once

#include <vector>
#include <sdsl/coder_elias_delta.hpp>

using namespace std;

class Delta_Vector{

public:

    sdsl::int_vector<> data; // Elias-delta-encoded diffs data[i] = v[i] - v[i-1], such that data[0] := v[0] + 1

    // The input values must be sorted and distinct
    Delta_Vector(const vector<int64_t>& sorted_values){
        sdsl::int_vector<> diffs(sorted_values.size());
        diffs[0] = sorted_values[0] + 1; // Elias-delta does not have a code for zero, so we add one
        for(int64_t i = 1; i < sorted_values.size(); i++)
            diffs[i] = sorted_values[i] - sorted_values[i-1];
        sdsl::coder::elias_delta::encode(diffs, data);
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

};