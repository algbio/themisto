#pragma once

#include "globals.hh"
#include "sbwt/globals.hh"
#include "sbwt/buffered_streams.hh"
#include "sbwt/EM_sort/EM_sort.hh"
#include "sbwt/EM_sort/bit_level_stuff.hh"
#include "sdsl/bit_vectors.hpp"


// Should be constructed with Sparse_Uint_Array_Builder
class Sparse_Uint_Array{
    private:

    Sparse_Uint_Array(const Sparse_Uint_Array& temp_obj) = delete; // No copying because pointers inside
    Sparse_Uint_Array& operator=(const Sparse_Uint_Array& temp_obj) = delete;  // No copying because pointers inside

    sdsl::bit_vector marks; // Marks which cells have a value
    sdsl::rank_support_v5<> marks_rs; // Rank support for marks
    sdsl::int_vector<> values; // Values at those cells that are marked
    uint64_t max_value;

    public:

    Sparse_Uint_Array(){}

    Sparse_Uint_Array(const sdsl::bit_vector& marks, sdsl::int_vector<>& values, uint64_t max_value) :
        marks(marks), values(values), max_value(max_value){
            sdsl::util::init_support(marks_rs, &marks);
        }


    // Return -1 if not in the array
    int64_t get(uint64_t idx) const{
        if(idx >= marks.size()) throw std::runtime_error("Access out of bounds at Sparse_Uint_Array");
        if(marks[idx] == 0) return -1; // Not in array
        int64_t pos = marks_rs.rank(idx);
        return values[pos];
    }

    bool has_index(uint64_t idx) const{
        return marks[idx];
    }

    int64_t serialize(ostream& os) const{
        int64_t n_bytes_written = 0;
        n_bytes_written += marks.serialize(os);
        n_bytes_written += marks_rs.serialize(os);
        n_bytes_written += values.serialize(os);

        os.write((char*)&max_value, sizeof(max_value));
        n_bytes_written += sizeof(max_value);

        return n_bytes_written;
    }

    void load(istream& is){
        marks.load(is);
        marks_rs.load(is);
        values.load(is);
        is.read((char*)&max_value, sizeof(max_value));

        marks_rs.set_vector(&marks);
    }
};

class Sparse_Uint_Array_Builder{

    private:

    string temp_filename;
    sbwt::Buffered_ofstream<> out_stream;
    sdsl::bit_vector marks; // Marks which indices have been set
    uint64_t array_length;
    uint64_t max_value;
    uint64_t n_values;
    uint64_t ram_bytes;
    uint64_t n_threads;

    string EM_sort_big_endian_LL_pairs(string infile, uint64_t ram_bytes, uint64_t key, uint64_t n_threads){
        assert(key == 0 || key == 1);
        
        auto cmp = [&](const char* A, const char* B) -> bool{
            uint64_t x_1,y_1,x_2,y_2; // compare (x_1, y_1) and (x_2,y_2)
            x_1 = sbwt::parse_big_endian_LL(A + 0);
            y_1 = sbwt::parse_big_endian_LL(A + 8);
            x_2 = sbwt::parse_big_endian_LL(B + 0);
            y_2 = sbwt::parse_big_endian_LL(B + 8);
            if(key == 0){
                return make_pair(x_1, y_1) < make_pair(x_2, y_2);
            } else{
                return make_pair(y_1, x_1) < make_pair(y_2, x_2);
            }
        };

        string outfile = sbwt::get_temp_file_manager().create_filename();
        sbwt::EM_sort_constant_binary(infile, outfile, cmp, ram_bytes, 8+8, n_threads);
        return outfile;
    }


    public:

    Sparse_Uint_Array_Builder(uint64_t array_length, uint64_t ram_bytes, uint64_t n_threads) : array_length(array_length), max_value(0), n_values(0), ram_bytes(ram_bytes), n_threads(n_threads) {
        temp_filename = sbwt::get_temp_file_manager().create_filename("");
        out_stream.open(temp_filename, ios::binary);
        marks = sdsl::bit_vector(array_length, 0);
    }

    // If multiple values are added to the same index, the *smallest* value to that index is kept
    void add(uint64_t index, uint64_t value){
        assert(value <= max_value);
        marks[index] = 1;
        write_big_endian_LL(out_stream, index);
        write_big_endian_LL(out_stream, value);
        max_value = max(max_value, value);
        n_values++;
    }    

    Sparse_Uint_Array finish(){

        out_stream.close();

        string sorted_out = EM_sort_big_endian_LL_pairs(temp_filename, ram_bytes, 0, n_threads);
        sbwt::get_temp_file_manager().delete_file(temp_filename);
        sbwt::Buffered_ifstream<> sorted_in(sorted_out);
        vector<char> buffer(8+8);
        uint64_t rank = 0;
        int64_t prev_index = -1;
        sdsl::int_vector<> values(n_values, 0, ceil(log2(max_value)));
        while(true){
            sorted_in.read(buffer.data(), 8+8);
            if(sorted_in.eof()) break;
            uint64_t index = sbwt::parse_big_endian_LL(buffer.data());
            uint64_t value = sbwt::parse_big_endian_LL(buffer.data() + 8);

            if(index != prev_index) values[rank++] = value;
            // If there are multiple values with the same index, we ignore all but the first one
            // This is why the comment on add(...) says that the smallest value is kept

            prev_index = index;
        }
        sbwt::get_temp_file_manager().delete_file(sorted_out);

        return Sparse_Uint_Array(marks, values, max_value);
    }

};
