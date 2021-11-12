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
#include "globals.hh"
#include "Block.hh"
#include "ParallelBoundedQueue.hh"
#include "generic_EM_classes.hh"

/*
Design:
  
There are two operating modes: EM_LINES and EM_BINARY. The mode determines
how records are represented in disk and memory.

EM_LINES:
  - The input is a file with one line per record. user-provided comparison function 
    takes as a parameter two NULL-terminated C-style strings which are lines in the
    input file.
  - Internal workings: The records are read from disk into an object of class Block. The 
    data of the block is stored in one large char-array. All lines are concatenated, 
    but newlines are replaced with NULL-terminator so that strcmp can be used easily
    for doing the comparisons if wanted. A vector<LL> 'starts' stores the 
    start positions of each record. Sorting is done by permuting 'starts' and keeping the 
    data in place. When the data is written back to diks, the null-terminators are again 
    replaced by newline-characters, so the resulting file is a permutation of the lines 
    of the original file.

EM_BINARY:
  - Records are in binary such that the first 8 bytes are the big-endian representation
    of a 64-bit integer x. Then follows x-8 bytes of data. The comparison function now takes
    as a parameter the pointers to the starts of the x byte block of data of the records.
    When passing around records, the size x is always at the start of the data array.
  - Internal workings: The records are read from disk into an object of class Block, which
    has a char-array which will store the concatenation of all records. A vector<LL> 'starts'
    stores the start position of each record. Sorting is done by permuting 'starts' and keeping
    the data in place. The data is written back to disk according to the permutation of starts.
*/

using namespace std;

const int EM_LINES = 0;
const int EM_VARIABLE_BINARY = 1;
const int EM_CONSTANT_BINARY = 2;

// Interprets the strings as integers (no leading zeros allowed) and returns:
//     -1 if x < y
//      0 if x = y
//      1 if x > y
int compare_as_numbers(const char* x, const char* y){
    LL nx = strlen(x);
    LL ny = strlen(y);
    if(nx < ny) return -1;
    if(nx > ny) return 1;
    return strcmp(x,y);
}

bool memcmp_variable_binary_records(const char* x, const char* y){
    LL nx = parse_big_endian_LL(x);
    LL ny = parse_big_endian_LL(y);
    LL c = memcmp(x + 8, y + 8, min(nx-8,ny-8));
    if(c < 0){
        return true;
    }
    else if(c > 0){
        return false;
    }
    else { // c == 0
        return nx < ny;
    }
}

void copy_file(string infile, string outfile, LL buf_size){
    throwing_ifstream in(infile, ios::binary);
    throwing_ofstream out(outfile, ios::binary);

    vector<char> buf(buf_size);
    while(true){
        in.read(buf.data(), buf_size);
        LL bytes_read = in.gcount();
        if(bytes_read > 0){
            out.write(buf.data(), bytes_read);
        } else break;
    }
}

void merge_files_generic(const std::function<bool(const char* x, const char* y)>& cmp, LL& merge_count, Generic_Record_Reader* reader, Generic_Record_Writer* writer){

    write_log("Doing merge number " + to_string(merge_count));

    vector<char*> input_buffers;
    vector<LL> input_buffer_sizes;

    for(int64_t i = 0; i < reader->get_num_files(); i++){
        LL buf_size = 1024;
        input_buffers.push_back((char*)malloc(buf_size)); // Freed at the end of this function
        input_buffer_sizes.push_back(buf_size);
    }

    auto cmp_wrap = [&](pair<char*, int64_t> x, pair<char*, int64_t> y) { 
        return cmp(x.first, y.first);
    };

    multiset<pair<char*, int64_t>, decltype(cmp_wrap)> Q(cmp_wrap); // Priority queue: (record, file index).
    // Must be a multiset because a regular set will not add an element if it is equal
    // to another according to the comparison function.

    // Initialize priority queue
    string line;
    for(int64_t i = 0; i < reader->get_num_files(); i++){
        reader->read_record(i, &input_buffers[i], &input_buffer_sizes[i]);
        Q.insert({input_buffers[i], i});
    }

    // Do the merge
    while(!Q.empty()){

        char* record; LL stream_idx;
        std::tie(record, stream_idx) = *(Q.begin());
        Q.erase(Q.begin()); // pop

        // Write the current data
        writer->write(record);

        // Read next value from the file
        if(reader->read_record(stream_idx, &input_buffers[stream_idx], &input_buffer_sizes[stream_idx]))
            Q.insert({input_buffers[stream_idx], stream_idx});
    }

    writer->close_file();

    for(int64_t i = 0; i < reader->get_num_files(); i++){
        free(input_buffers[i]);
    }

    merge_count++;

}

void EM_sort_generic(string infile, string outfile, const std::function<bool(const char* x, const char* y)>& cmp, LL RAM_bytes, LL k, Generic_Block_Producer* producer, vector<Generic_Block_Consumer*> consumers, Generic_Record_Reader* reader, Generic_Record_Writer* writer){
    assert(k >= 2);

    // Number of blocks in the memory at once:
    // - 1 per consumer thread in processing
    // - 1 in the queue
    // - 1 with the producer loading (there is only one producer)
    // So if block size is B, we have (n_threads + 2)*B blocks in memory at a time
    // So we have the equation (n_threads + 2)*B = RAM_bytes. Solve for B:

    LL B = RAM_bytes / (consumers.size() + 2);

    vector<string> block_files;
    ParallelBoundedQueue<Generic_Block*> Q(1); // 1 byte = basically only one block can be in the queue at a time
    vector<std::thread> threads;
    
    // Create producer
    threads.push_back(std::thread([&Q,&infile,&B,&producer](){
        producer->run(Q, B);
    }));

    // Create consumers
    for(int64_t i = 0; i < consumers.size(); i++){
        threads.push_back(std::thread([i, &Q, &block_files, &cmp, &consumers](){
            write_log("Starting thread " + to_string(i));
            consumers[i]->run(Q,cmp);
            write_log("Thread " + to_string(i) + ": done");
        }));
    }

    for(std::thread& t : threads) t.join();

    for(auto& consumer : consumers){
        for(string filename : consumer->get_outfilenames()){
            block_files.push_back(filename);
        }
    }

    // Merge blocks
    LL merge_count = 0;
    vector<string> cur_round = block_files;
    while(cur_round.size() > 1){
        vector<string> next_round;
        for(LL i = 0; i < cur_round.size(); i += k){
            // Merge
            vector<string> to_merge(cur_round.begin() + i, cur_round.begin() + min(i + k, (LL)cur_round.size()));
            string round_file = temp_file_manager.get_temp_file_name("");
            writer->open_file(round_file);
            reader->open_files(to_merge);
            merge_files_generic(cmp, merge_count, reader, writer);
            next_round.push_back(round_file);
            writer->close_file();
            reader->close_files();

            // Clear files
            for(LL j = i; j < min(i+k, (LL)cur_round.size()); j++){
                temp_file_manager.delete_file(cur_round[j].c_str());
            }
        }
        cur_round = next_round;
    }

    // Copy final merge file to outfile
    if(cur_round.size() == 0) // Function was called with empty input file
        copy_file(infile, outfile, 1024*1024);
    else{
        assert(cur_round.size() == 1);
        copy_file(cur_round[0], outfile, 1024*1024);
        temp_file_manager.delete_file(cur_round[0].c_str());
    }

}

// Constant size records of record_size bytes each
void EM_sort_constant_binary(string infile, string outfile, const std::function<bool(const char* x, const char* y)>& cmp, LL RAM_bytes, LL k, LL record_size, LL n_threads){
    assert(k >= 2);

    Generic_Block_Producer* producer = new Constant_Block_Producer(infile, record_size);
    vector<Generic_Block_Consumer*> consumers;
    for(LL i = 0; i < n_threads; i++)
        consumers.push_back(new Block_Consumer(i));
    Generic_Record_Reader* reader = new Constant_Record_Reader(record_size);
    Generic_Record_Writer* writer = new Constant_Record_Writer(record_size);

    EM_sort_generic(infile, outfile, cmp, RAM_bytes, k, producer, consumers, reader, writer);

    delete producer;
    for(Generic_Block_Consumer* C : consumers) delete C;
    delete reader;
    delete writer;

}

// Binary format of record: first 8 bytes give the length of the record, then comes the record
// Line format of record: A char sequence terminated by \n\0. For the cmp function the terminator is just \0.
// k = k-way merge parameter
// mode = EM_LINES or mode = EM_VARIABLE_BINARY or mode = EM_CONSTANT_BINARY
void EM_sort(string infile, string outfile, const std::function<bool(const char* x, const char* y)>& cmp, LL RAM_bytes, LL k, LL n_threads, LL mode){ // todo: rename to EM_sort_variable

    assert((mode == EM_LINES || mode == EM_VARIABLE_BINARY) && k >= 2);

    Generic_Block_Producer* producer; 
    if(mode == EM_LINES)
        producer = new Line_Block_Producer(infile);
    if(mode == EM_VARIABLE_BINARY)
        producer = new Variable_Block_Producer(infile);

    vector<Generic_Block_Consumer*> consumers;

    for(LL i = 0; i < n_threads; i++)
        consumers.push_back(new Block_Consumer(i));
    
    Generic_Record_Reader* reader;
    if(mode == EM_LINES)
        reader = new Line_Record_Reader();
    if(mode == EM_VARIABLE_BINARY)
        reader = new Variable_Record_Reader();

    Generic_Record_Writer* writer;
    if(mode == EM_LINES)
        writer = new Line_Record_Writer();
    if(mode == EM_VARIABLE_BINARY)
        writer = new Variable_Record_Writer();

    EM_sort_generic(infile, outfile, cmp, RAM_bytes, k, producer, consumers, reader, writer);

    delete producer;
    for(Generic_Block_Consumer* C : consumers) delete C;
    delete reader;
    delete writer;

}
