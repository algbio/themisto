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
#include <cassert>
#include <functional>
#include "globals.hh"
#include "Block.hh"
#include "ParallelBoundedQueue.hh"


class Generic_Block_Producer{
public:
    virtual void run(ParallelBoundedQueue<Generic_Block*>& Q, LL block_size) = 0;
    virtual ~Generic_Block_Producer(){}
};

class Generic_Block_Consumer{
public:
    virtual void run(ParallelBoundedQueue<Generic_Block*>& Q, const std::function<bool(const char* x, const char* y)>& cmp) = 0;
    virtual vector<string> get_outfilenames() = 0;
    virtual ~Generic_Block_Consumer(){}
};

class Generic_Record_Reader{
public:

    virtual void open_files(vector<string> filenames) = 0;
    virtual void close_files() = 0;
    virtual LL get_num_files() = 0;
    virtual bool read_record(LL input_index, char** buffer, LL* buffer_size) = 0; // reallocs. Returns true iff read succeeded
    virtual ~Generic_Record_Reader(){}

};

class Generic_Record_Writer{
public:

    virtual void write(char* record) = 0;
    virtual void open_file(string filename) = 0;
    virtual void close_file() = 0;
    virtual ~Generic_Record_Writer(){}

};

class Block_Consumer : public Generic_Block_Consumer{ // Works for any type of block
public:

    vector<string> filenames;
    LL thread_id;

    Block_Consumer(LL thread_id) : thread_id(thread_id) {}

    virtual void run(ParallelBoundedQueue<Generic_Block*>& Q, const std::function<bool(const char* x, const char* y)>& cmp){
        while(true){
            Generic_Block* block  = Q.pop();
            if(block == nullptr){
                Q.push(nullptr, 0);
                break; // No more work available
            }
            write_log("Thread " + to_string(thread_id) + ": Starting to sort a block.");
            block->sort(cmp);
            write_log("Thread " + to_string(thread_id) + ": Finished sorting a block.");
            filenames.push_back(create_temp_filename());
            block->write_to_file(filenames.back());
            delete block; // Initially allocated by a producer
        }
    }

    virtual vector<string> get_outfilenames(){
        return filenames;
    }
};

//
// Constant record classes below
//

class Constant_Block_Producer : public Generic_Block_Producer{
    public:

    throwing_ifstream in;
    LL record_size;

    Constant_Block_Producer(string infile, LL record_size) : in(infile, ios::binary), record_size(record_size) {}

    virtual void run(ParallelBoundedQueue<Generic_Block*>& Q, LL block_size){
        while(true){
            Constant_binary_block* block = get_next_constant_binary_block(in,block_size,record_size); // Freed by a consumer
            if(block->starts.size() == 0){
                Q.push(nullptr, 0);
                break; // No more work available
            }
            Q.push((Generic_Block*)block, block_size);
        }
    }
};

class Constant_Record_Reader : public Generic_Record_Reader{
public:

    vector<throwing_ifstream> inputs;
    LL record_size;

    Constant_Record_Reader(LL record_size) : record_size(record_size){}

    virtual void open_files(vector<string> filenames){
        inputs.clear();
        inputs.resize(filenames.size());
        for(LL i = 0; i < filenames.size(); i++){
            inputs[i].open(filenames[i], ios::binary);
        }
    }

    virtual void close_files(){
        for(throwing_ifstream& in : inputs) in.close();
    }

    virtual LL get_num_files(){
        return inputs.size();
    }

    virtual bool read_record(LL input_index, char** buffer, LL* buffer_size){
        if(*buffer_size < record_size){
            *buffer = (char*)realloc(*buffer, record_size);
            *buffer_size = record_size;
        }
        return inputs[input_index].read(*buffer, record_size);
    }

};

class Constant_Record_Writer : public Generic_Record_Writer{
public:
    throwing_ofstream out;
    LL record_size;

    Constant_Record_Writer(LL record_size) : record_size(record_size){}

    virtual void open_file(string filename){
        out.open(filename, ios::binary);
    }

    virtual void close_file(){
        out.close();
    }

    virtual void write(char* record){
        out.write(record, record_size);
    }

};

//
// Variable binary record classes below
//

class Variable_Block_Producer : public Generic_Block_Producer{
    public:

    throwing_ifstream in;

    Variable_Block_Producer(string infile) : in(infile, ios::binary) {}

    virtual void run(ParallelBoundedQueue<Generic_Block*>& Q, LL block_size){
        while(true){
            Variable_binary_block* block = get_next_variable_binary_block(in,block_size); // Freed by a consumer
            if(block->starts.size() == 0){
                Q.push(nullptr, 0);
                break; // No more work available
            }
            Q.push((Generic_Block*)block, block_size);
        }
    }
};

class Variable_Record_Reader : public Generic_Record_Reader{
public:

    vector<throwing_ifstream> inputs;

    Variable_Record_Reader() {}

    virtual void open_files(vector<string> filenames){
        inputs.clear();
        inputs.resize(filenames.size());
        for(LL i = 0; i < filenames.size(); i++){
            inputs[i].open(filenames[i], ios::binary);
        }
    }

    virtual void close_files(){
        for(throwing_ifstream& in : inputs) in.close();
    }

    virtual LL get_num_files(){
        return inputs.size();
    }

    virtual bool read_record(LL input_index, char** buffer, LL* buffer_size){
        return read_variable_binary_record(inputs[input_index], buffer, buffer_size);
    }

};

class Variable_Record_Writer : public Generic_Record_Writer{
public:
    throwing_ofstream out;

    Variable_Record_Writer() {}

    virtual void open_file(string filename){
        out.open(filename, ios::binary);
    }

    virtual void close_file(){
        out.close();
    }

    virtual void write(char* record){
        LL rec_len = parse_big_endian_LL(record);
        out.write(record, rec_len);
    }

};


//
// Line record classes below
//


class Line_Block_Producer : public Generic_Block_Producer{
public:

    throwing_ifstream in;

    Line_Block_Producer(string infile) : in(infile) {}

    virtual void run(ParallelBoundedQueue<Generic_Block*>& Q, LL block_size){
        while(true){
            Line_block* block = get_next_line_block(in,block_size); // Freed by a consumer
            if(block->starts.size() == 0){
                Q.push(nullptr, 0);
                break; // No more work available
            }
            Q.push(block, block_size);
        }
    }
};

class Line_Record_Reader : public Generic_Record_Reader{
public:

    vector<throwing_ifstream> inputs;
    string line; // Reusable buffer

    Line_Record_Reader() {}

    virtual void open_files(vector<string> filenames){
        inputs.clear();
        inputs.resize(filenames.size());
        for(LL i = 0; i < filenames.size(); i++){
            inputs[i].open(filenames[i]);
        }
    }

    virtual void close_files(){
        for(throwing_ifstream& in : inputs) in.close();
    }

    virtual LL get_num_files(){
        return inputs.size();
    }

    virtual bool read_record(LL input_index, char** buffer, LL* buffer_size){
        if(inputs[input_index].getline(line)){
            while(*buffer_size < line.size() + 1){ // +1: null terminator
                *buffer = (char*)realloc(*buffer, (*buffer_size)*2);
                *buffer_size *= 2;
            }
            memcpy(*buffer, line.c_str(), line.size() + 1); // '\n' is now replaced with \0
            return true;
        }
        return false;
    }

};

class Line_Record_Writer : public Generic_Record_Writer{
public:
    throwing_ofstream out;

    Line_Record_Writer() {}

    virtual void open_file(string filename){
        out.open(filename);
    }

    virtual void close_file(){
        out.close();
    }

    virtual void write(char* record){
        LL len = 0;
        while(*(record+len) != 0) len++;
        out.write(record, len);
        out.write("\n", 1);
    }

};

