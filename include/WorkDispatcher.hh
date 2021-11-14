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

struct Read_Batch{
    // This is better than vector<string> because this does not have memory fragmentation

    int64_t read_id_of_first; // read id of first read in the batch
    string data; // concatenation of reads
    vector<bool> read_starts; // marks start of reads
};

class Batch_reader{
// An iterator reading one read starting from the given start position

public:

    Read_Batch* batch;
    int64_t pos; // Current position

    Batch_reader(Read_Batch* batch, int64_t start) : batch(batch), pos(start) {}

    // Returns (pointer to start, length). If done, return (NULL, 0)
    pair<const char*, LL> get_next_string(){
        if(pos >= batch->data.size()) return {NULL,0}; // End of batch
        LL len = 1;
        LL start = pos;
        while(true){
            pos++;
            if(pos >= batch->data.size() || batch->read_starts[pos] == 1) break;
            len++;
        }
        return {batch->data.c_str() + start, len};
    }
};

class DispatcherConsumerCallback{
public:
    virtual void callback(const char* S, LL S_size, int64_t string_id) = 0;
    virtual void finish() = 0;
    virtual ~DispatcherConsumerCallback() {} 
};

void dispatcher_consumer(ParallelBoundedQueue<Read_Batch>& Q, DispatcherConsumerCallback* cb, LL thread_id);

// Will run characters through fix_char, which at the moment of writing this comment
// upper-cases the character and further the result is not A, C, G or T, changes it to A.
void dispatcher_producer(ParallelBoundedQueue<Read_Batch>& Q, Sequence_Reader& sr, int64_t buffer_size);

void run_dispatcher(vector<DispatcherConsumerCallback*>& callbacks, Sequence_Reader& sr, LL buffer_size);