#pragma once
#include <string>
#include <vector>
#include <memory>
#include <array>

#include "globals.hh"
#include "zstr.hpp"
#include "old_buffered_streams.hh"
#include "sbwt/EM_sort/ParallelBoundedQueue.hh"

using namespace std;
using namespace sbwt;

class ParallelBaseWriter{

public:

virtual void write(const string& data) = 0;
virtual void write(const char* data, int64_t data_length) = 0;
virtual void flush() = 0;
virtual ~ParallelBaseWriter(){}

};

class ParallelOutputWriter : public ParallelBaseWriter{
    public:

    string outfile;
    old_themisto_code::Buffered_ofstream outstream;
    std::mutex mutex;

    ParallelOutputWriter(string outfile) : outfile(outfile){
        outstream.open(outfile);
    }


    ParallelOutputWriter(ostream &ostream) : outstream(ostream) {}


    virtual void write(const string& data){
        write(data.data(), data.size());
    }

    virtual void write(const char* data, int64_t data_length){
        std::lock_guard<std::mutex> lg(mutex);
        outstream.write(data, data_length);        
    }

    virtual void flush(){
        outstream.flush();
    }
    
};

// TODO: is not buffered
class ParallelGzipWriter : public  ParallelBaseWriter{
    public:

    string outfile;
    unique_ptr<zstr::ostream> gzip_outstream;
    std::mutex mutex;

    ParallelGzipWriter(string outfile) : outfile(outfile){
        check_writable(outfile);
	gzip_outstream = unique_ptr<zstr::ostream>(reinterpret_cast<zstr::ostream*>(new zstr::ofstream(outfile)));
    }

    ParallelGzipWriter(ostream &ostream) {
	gzip_outstream = unique_ptr<zstr::ostream>(new zstr::ostream(ostream.rdbuf()));
    }

    virtual void write(const string& data){
        write(data.data(), data.size());
    }

    virtual void write(const char* data, int64_t data_length){
        std::lock_guard<std::mutex> lg(mutex);
        gzip_outstream.get()->write(data, data_length);
    }

    virtual void flush(){
        gzip_outstream->flush();
    }

    virtual ~ParallelGzipWriter(){
        gzip_outstream.reset();
    }
};

// Todo: this class is no longer needed because ParallelBaseWriter now has a write-function with c-string input
class ParallelBinaryOutputWriter{
    public:

    string outfile;
    old_themisto_code::Buffered_ofstream outstream;
    std::mutex mutex;

    ParallelBinaryOutputWriter(string outfile) : outfile(outfile, ios::binary){
        outstream.open(outfile);
    }

    void write(const char* data, int64_t n_bytes){
        std::lock_guard<std::mutex> lg(mutex);
        outstream.write(data, n_bytes);
    }

    void flush(){
        outstream.flush();
    }
};


struct ReadBatch{
    public:
    uint64_t firstReadID; // read id of first read in the batch
    vector<char> data; // concatenation of reads
    vector<uint64_t> readStarts; // indicates starting positions of reads in data array, last item is a dummy
    vector<std::array<uint8_t, 8>> metadata; // For each read, 8 bytes of metadata

    int64_t byte_size() const{
        return sizeof(firstReadID) * 1 + 
               sizeof(char) * data.size() + 
               sizeof(uint64_t) * readStarts.size() + 
               sizeof(uint8_t) * 8 * metadata.size();
    }
};

class ReadBatchIterator{
// An iterator reading one read starting from the given start position

public:

    ReadBatch* batch;
    uint64_t pos; // Index in the readStarts array of the batch
    std::array<uint8_t, 8> dummy_metadata;

    ReadBatchIterator(ReadBatch* batch, int64_t start) : batch(batch), pos(start) {}

    // Returns (c-string, length, metadata). If done, return (NULL, 0, NULL)
    std::tuple<const char*, uint64_t, std::array<uint8_t, 8>> getNextRead(){
        if(pos >= (batch->readStarts.size()-1)) return {NULL, 0, dummy_metadata}; // End of batch
        uint64_t len = batch->readStarts[pos+1] - batch->readStarts[pos];
        uint64_t start = batch->readStarts[pos];
        std::array<uint8_t, 8> metadata = batch->metadata[pos];
        pos++; // Move to next read ready for next call to this function
        return {(batch->data.data()) + start, len, metadata};
    }
};

class DispatcherConsumerCallback{
public:

    // Takes a string S (not necessarily null-terminated), the length of S, the read_id of S and 8 bytes of metadata.
    // It's 8 bytes so that in the future if we want to pass in arbitrary data, it can fit a 8-byte pointer to memory.
    // The main use case currently is to pass a 64-bit color id. Why not just pass a void pointer to the color id? 
    // Because then the memory management gets unnecessarily complicated and inefficient in the 64-bit integer use case.
    virtual void callback(const char* S, int64_t S_size, int64_t read_id, std::array<uint8_t, 8> metadata) = 0;
    virtual void finish() = 0;
    virtual ~DispatcherConsumerCallback() {} 
};

class Metadata_Stream{

    public:

    virtual std::array<uint8_t, 8> next() = 0; // Undefined behavior if called after all the metadata has been read

};


void dispatcher_consumer(ParallelBoundedQueue<ReadBatch*>& Q, DispatcherConsumerCallback* cb, int64_t thread_id);

// Will run characters through fix_char, which at the moment of writing this comment
// upper-cases the character and further the result is not A, C, G or T, changes it to A.
template<typename sequence_reader_t>
void dispatcher_producer(ParallelBoundedQueue<ReadBatch*>& Q, sequence_reader_t& sr, Metadata_Stream* metadata_stream, int64_t batch_size){
    // Push work in batches of approximately buffer_size base pairs

    int64_t read_id = 0;
    ReadBatch* batch = new ReadBatch(); // Deleted by consumer
    std::array<uint8_t, 8> dummy_metadata;

    auto push_batch = [&](){
        batch->readStarts.push_back(batch->data.size()); // Append the end sentinel
        batch->metadata.push_back(dummy_metadata); // Append the end sentinel
        Q.push(batch, batch->byte_size());
        batch = new ReadBatch(); // Clear
    };

    while(true){
        // Create a new read batch and push it to the work queue

        int64_t len = sr.get_next_read_to_buffer();
        if(len == 0){
            // All reads read. Push the current batch if it's non-empty and quit
            if(batch->data.size() > 0) push_batch();
            break; // quit
        } else{
            // Add read to batch
            if(batch->data.size() == 0) batch->firstReadID = read_id;
            batch->readStarts.push_back(batch->data.size());
            batch->metadata.push_back(metadata_stream->next());
            for(int64_t i = 0; i < len; i++) batch->data.push_back(sr.read_buf[i]);
            if(batch->data.size() >= batch_size) push_batch();
            read_id++;
        }
    }
    
    batch->readStarts.push_back(batch->data.size()); // Append the end sentinel
    batch->metadata.push_back(dummy_metadata); // Append the end sentinel
    Q.push(batch,0); // Empty batch in the end signifies end of the queue
}

template<typename sequence_reader_t>
void run_dispatcher(vector<DispatcherConsumerCallback*>& callbacks, sequence_reader_t& sr, Metadata_Stream* metadata_stream, int64_t buffer_size){
    vector<std::thread> threads;
    ParallelBoundedQueue<ReadBatch*> Q(buffer_size);

    // Create consumers
    for(int64_t i = 0; i < callbacks.size(); i++){
        threads.push_back(std::thread(dispatcher_consumer,std::ref(Q),callbacks[i], i));
    }

    // Create producer
    threads.push_back(std::thread(dispatcher_producer<sequence_reader_t>,std::ref(Q),std::ref(sr),metadata_stream,buffer_size));

    for(std::thread& t : threads) t.join();

    for(int64_t i = 0; i < callbacks.size(); i++){
        callbacks[i]->finish();
    }
}
