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
#include "EM_sort.hh"

using namespace std;

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

template <typename record_reader_t, typename record_writer_t>
void merge_files_generic(const std::function<bool(const char* x, const char* y)>& cmp, LL& merge_count, record_reader_t& reader, record_writer_t& writer){

    write_log("Doing merge number " + to_string(merge_count));

    vector<char*> input_buffers;
    vector<LL> input_buffer_sizes;

    for(int64_t i = 0; i < reader.get_num_files(); i++){
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
    for(int64_t i = 0; i < reader.get_num_files(); i++){
        reader.read_record(i, &input_buffers[i], &input_buffer_sizes[i]);
        Q.insert({input_buffers[i], i});
    }

    // Do the merge
    while(!Q.empty()){

        char* record; LL stream_idx;
        std::tie(record, stream_idx) = *(Q.begin());
        Q.erase(Q.begin()); // pop

        // Write the current data
        writer.write(record);

        // Read next value from the file
        if(reader.read_record(stream_idx, &input_buffers[stream_idx], &input_buffer_sizes[stream_idx]))
            Q.insert({input_buffers[stream_idx], stream_idx});
    }

    writer.close_file();

    for(int64_t i = 0; i < reader.get_num_files(); i++){
        free(input_buffers[i]);
    }

    merge_count++;

}

template <typename record_reader_t, typename record_writer_t>
void EM_sort_generic(string infile, string outfile, const std::function<bool(const char* x, const char* y)>& cmp, LL RAM_bytes, Generic_Block_Producer* producer, vector<Generic_Block_Consumer*> consumers, record_reader_t& reader, record_writer_t& writer){

    LL max_files = 512;

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
        for(LL i = 0; i < cur_round.size(); i += max_files){
            // Merge
            vector<string> to_merge(cur_round.begin() + i, cur_round.begin() + min(i + max_files, (LL)cur_round.size()));
            string round_file = get_temp_file_manager().create_filename();
            writer.open_file(round_file);
            reader.open_files(to_merge);
            merge_files_generic(cmp, merge_count, reader, writer);
            next_round.push_back(round_file);
            writer.close_file();
            reader.close_files();

            // Clear files
            for(LL j = i; j < min(i+max_files, (LL)cur_round.size()); j++){
                get_temp_file_manager().delete_file(cur_round[j].c_str());
            }
        }
        cur_round = next_round;
    }

    // Move final merge file to outfile
    
    if(cur_round.size() == 0) // Function was called with empty input file
        std::filesystem::rename(infile, outfile);
    else{
        assert(cur_round.size() == 1);
        std::filesystem::rename(cur_round[0], outfile);
        get_temp_file_manager().delete_file(cur_round[0].c_str());
    }

}

// Constant size records of record_size bytes each
void EM_sort_constant_binary(string infile, string outfile, const std::function<bool(const char* x, const char* y)>& cmp, LL RAM_bytes, LL record_size, LL n_threads){

    Generic_Block_Producer* producer = new Constant_Block_Producer(infile, record_size);
    vector<Generic_Block_Consumer*> consumers;
    for(LL i = 0; i < n_threads; i++)
        consumers.push_back(new Block_Consumer(i));
    Constant_Record_Reader reader(record_size);
    Constant_Record_Writer writer(record_size);

    EM_sort_generic(infile, outfile, cmp, RAM_bytes, producer, consumers, reader, writer);

    delete producer;
    for(Generic_Block_Consumer* C : consumers) delete C;

}

void EM_sort_variable_length_records(string infile, string outfile, const std::function<bool(const char* x, const char* y)>& cmp, LL RAM_bytes, LL n_threads){

    Generic_Block_Producer* producer = new Variable_Block_Producer(infile);

    vector<Generic_Block_Consumer*> consumers;

    for(LL i = 0; i < n_threads; i++)
        consumers.push_back(new Block_Consumer(i));
    
    Variable_Record_Reader reader;
    Variable_Record_Writer writer;

    EM_sort_generic(infile, outfile, cmp, RAM_bytes, producer, consumers, reader, writer);

    delete producer;
    for(Generic_Block_Consumer* C : consumers) delete C;;

}
