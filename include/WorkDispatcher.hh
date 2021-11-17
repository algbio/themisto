#pragma once
#include <string>
#include <vector>
#include "ParallelBoundedQueue.hh"
#include "input_reading.hh"
#include "globals.hh"
#include "zstr.hpp"

using namespace std;

class ParallelBaseWriter{

public:

virtual void write(const string& result) = 0;
virtual void flush() = 0;
virtual ~ParallelBaseWriter(){}

};

class ParallelOutputWriter : public ParallelBaseWriter{
    public:

    string outfile;
    throwing_ofstream outstream;
    std::mutex mutex;

    ParallelOutputWriter(string outfile) : outfile(outfile){
        outstream.open(outfile);
    }

    virtual void write(const string& result){
        std::lock_guard<std::mutex> lg(mutex);
        outstream << result;
    }

    virtual void flush(){
        outstream.flush();
    }
    
};

class ParallelGzipWriter : public  ParallelBaseWriter{
    public:

    string outfile;
    zstr::ofstream* gzip_outstream;
    std::mutex mutex;

    ParallelGzipWriter(string outfile) : outfile(outfile){
        check_writable(outfile);
        gzip_outstream = new zstr::ofstream(outfile); 
    }

    virtual void write(const string& result){
        std::lock_guard<std::mutex> lg(mutex);
        *gzip_outstream << result;
    }

    virtual void flush(){
         gzip_outstream->flush();
    }

    virtual ~ParallelGzipWriter(){
        delete gzip_outstream;
    }
};

class ParallelBinaryOutputWriter{
    public:

    string outfile;
    throwing_ofstream outstream;
    std::mutex mutex;

    ParallelBinaryOutputWriter(string outfile) : outfile(outfile, ios::binary){
        outstream.open(outfile);
    }

    void write(const char* data, LL n_bytes){
        std::lock_guard<std::mutex> lg(mutex);
        outstream.write(data,n_bytes);
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
};

class ReadBatchIterator{
// An iterator reading one read starting from the given start position

public:

    ReadBatch* batch;
    uint64_t pos; // Index in the readStarts array of the batch

    ReadBatchIterator(ReadBatch* batch, int64_t start) : batch(batch), pos(start) {}

    // Returns (pointer to start, length). If done, return (NULL, 0)
    pair<const char*, uint64_t> getNextRead(){
        if(pos >= (batch->readStarts.size()-1)) return {NULL, 0}; // End of batch
        uint64_t len = batch->readStarts[pos+1] - batch->readStarts[pos];
        uint64_t start = batch->readStarts[pos];
        pos++; // Move to next read ready for next call to this function
        return {(batch->data.data()) + start, len};
    }
};

class DispatcherConsumerCallback{
public:
    virtual void callback(const char* S, LL S_size, int64_t string_id) = 0;
    virtual void finish() = 0;
    virtual ~DispatcherConsumerCallback() {} 
};

void dispatcher_consumer(ParallelBoundedQueue<ReadBatch*>& Q, DispatcherConsumerCallback* cb, LL thread_id);

// Will run characters through fix_char, which at the moment of writing this comment
// upper-cases the character and further the result is not A, C, G or T, changes it to A.
void dispatcher_producer(ParallelBoundedQueue<ReadBatch*>& Q, Sequence_Reader_Buffered& sr, int64_t buffer_size);

void run_dispatcher(vector<DispatcherConsumerCallback*>& callbacks, Sequence_Reader_Buffered& sr, LL buffer_size);