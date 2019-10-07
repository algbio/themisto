#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <set>
#include <utility>
#include <tuple>
#include <cstring>
#include <cstdio>
#include "bit_level_stuff.hh"

using namespace std;


class Generic_Block{
    public:
    virtual void sort(const std::function<bool(const char* x, const char* y)>& cmp) = 0;
    virtual void write_to_file(string filename) = 0;
    virtual ~Generic_Block(){}
};

// Takes a pointer to a buffer. Read to the buffer the whole record including the 8-byte size value.
// The buffer might be realloc'd. Returns a bool whether read was successful.
// buffer_len must be greater than 0 (otherwise buffer doubling does not work)
bool read_variable_binary_record(ifstream& input, char** buffer, LL* buffer_len){
    assert(*buffer_len > 0);
    char rec_len_buf[8];
    input.read(rec_len_buf, 8); // Try to read the length of the record
    while(input.gcount() > 0){
        // Read was successful
        LL rec_len = parse_big_endian_LL(rec_len_buf);
        while(*buffer_len < rec_len){ // Make space in the buffer if needed
            *buffer = (char*)realloc(*buffer, *(buffer_len)*2);
            *buffer_len *= 2;
        }
        memcpy(*buffer, rec_len_buf, 8);
        input.read(*buffer + 8, rec_len - 8); // Read the payload
        return true;
    }
    return false;
}


// Block of records of variable length. The first 8 bytes of a record are a big-endian integer L 
// that tells size of the record. Then follow L-8 bytes which is the "payload" of the record.
class Variable_binary_block : public Generic_Block{
    public:

    char* data; // A buffer. Not necessarily all full of data. Contains concateated records
    LL data_len = 0;
    LL next_start = 0;
    
    vector<LL> starts;

    Variable_binary_block(){
        data = (char*) malloc(1 * sizeof(char));
        data_len = 1;
    }

    ~Variable_binary_block(){
        free(data);
    }

    void double_space(){
        data = (char*)realloc(data, 2 * data_len * sizeof(char));
        data_len *= 2;
    }

    virtual void sort(const std::function<bool(const char* x, const char* y)>& cmp){
        auto cmp_wrap = [&](LL x, LL y){
            return cmp(data+x,data+y);
        };
        std::sort(starts.begin(), starts.end(), cmp_wrap);
    }

    void add_record(const char* record){
        LL space_left = data_len - next_start;
        LL rec_len = parse_big_endian_LL(record);
        while(space_left < rec_len){
            double_space();
            space_left = data_len - next_start;
        }

        memcpy(data+next_start, record, rec_len);
        starts.push_back(next_start);
        next_start += rec_len;
    }

    virtual void write_to_file(string filename){
        ofstream out(filename, ios::binary);
        for(LL i = 0; i < starts.size(); i++){
            LL length = parse_big_endian_LL(data + starts[i]);
            out.write(data + starts[i], length);
        }
    }

    LL estimate_size_in_bytes(){
        return next_start + starts.size() * sizeof(LL);
    }
};

// RETURN VALUE MUST BE FREED BY CALLER
Variable_binary_block* get_next_variable_binary_block(ifstream& input, LL B){
    LL buffer_len = 1024; // MUST HAVE AT LEAST 8 BYTES
    char* buffer = (char*)malloc(buffer_len);
    Variable_binary_block* block = new Variable_binary_block();

    while(read_variable_binary_record(input, &buffer, &buffer_len)){
        block->add_record(buffer);
        if(block->estimate_size_in_bytes() > B) break;
    }

    free(buffer);
    return block;
}


// Records of exactly n bytes each
class Constant_binary_block : public Generic_Block{
public:

    public:

    char* data; // A buffer. Not necessarily all full of data. Contains concateated records
    LL data_len = 0;
    LL next_start = 0;
    LL record_size;
    vector<LL> starts; // We sort this vector. The data array is not touched in sorting.

    Constant_binary_block(LL record_size) : record_size(record_size) {
        data = (char*) malloc(1 * sizeof(char));
        data_len = 1;
    }

    ~Constant_binary_block(){
        free(data);
    }

    void double_space(){
        data = (char*)realloc(data, 2 * data_len * sizeof(char));
        data_len *= 2;
    }

    virtual void sort(const std::function<bool(const char* x, const char* y)>& cmp){
        auto cmp_wrap = [&](LL x, LL y){
            return cmp(data+x,data+y);
        };
        std::sort(starts.begin(), starts.end(), cmp_wrap);
    }

    void add_record(const char* record){
        LL space_left = data_len - next_start;
        while(space_left < record_size){
            double_space();
            space_left = data_len - next_start;
        }

        memcpy(data+next_start, record, record_size);
        starts.push_back(next_start);
        next_start += record_size;
    }

    virtual void write_to_file(string filename){
        ofstream out(filename, ios::binary);
        for(LL i = 0; i < starts.size(); i++){
            out.write(data + starts[i], record_size);
        }
    }

    LL estimate_size_in_bytes(){
        return next_start + starts.size() * sizeof(LL);
    }

};

// Reads up to B bytes into a new block
// THE RETURN VALUE MUST BE FREED BY THE CALLER
Constant_binary_block* get_next_constant_binary_block(ifstream& input, LL B, LL record_size){
    Constant_binary_block* block = new Constant_binary_block(record_size);
    char* buf = (char*)malloc(record_size);
    while(true){
        input.read(buf, record_size);
        if(input.gcount() == 0) break; // end of file
        block->add_record(buf);
        if(block->estimate_size_in_bytes() > B) break;
    }

    free(buf);
    return block;
}

// Records that are ascii-strings terminated with the \n character. In the internal
// data array the \n characters are replaced with null-terminators.
class Line_block : public Generic_Block{

    public:

    char* data; // A buffer. Not necessarily all full of data. Contains concateated records
    LL data_len = 0;
    LL next_start = 0;
    
    vector<LL> starts;

    Line_block(){
        data = (char*) malloc(1 * sizeof(char));
        data_len = 1;
    }

    ~Line_block(){
        free(data);
    }

    void double_space(){
        data = (char*)realloc(data, 2 * data_len * sizeof(char));
        data_len *= 2;
    }

    virtual void sort(const std::function<bool(const char* x, const char* y)>& cmp){
        auto cmp_wrap = [&](LL x, LL y){
            return cmp(data+x,data+y);
        };
        std::sort(starts.begin(), starts.end(), cmp_wrap);
    }

    void add_record(const string& line){
        LL space_left = data_len - next_start;
        while(space_left < line.size()+1){
            double_space();
            space_left = data_len - next_start;
        }
        starts.push_back(next_start);
        for(LL i = 0; i < line.size(); i++) {
            data[next_start + i] = line[i];
        }
        data[next_start + line.size()] = 0; // Null terminator
        next_start += line.size() + 1;
    }

    virtual void write_to_file(string filename){
        ofstream out(filename);
        for(LL i = 0; i < starts.size(); i++){
            string S(data + starts[i]);
            out << S << "\n";
        }
    }

    LL estimate_size_in_bytes(){
        return next_start + starts.size() * sizeof(LL);
    }

};

// Returns a heap object. MUST BE FREED BY THE CALLER
Line_block* get_next_line_block(ifstream& input, LL B){
    string line;
    Line_block* block = new Line_block();

    while(getline(input,line)){
        block->add_record(line);
        if(block->estimate_size_in_bytes() > B) break;
    }

    return block;
}
