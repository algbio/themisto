#pragma once
#include <string>
#include <vector>
#include "ParallelBoundedQueue.hh"
#include "input_reading.hh"
#include "globals.hh"

using namespace std;

class ParallelOutputWriter{
    public:

    string outfile;
    throwing_ofstream outstream;
    std::mutex mutex;

    ParallelOutputWriter(string outfile) : outfile(outfile){
        outstream.open(outfile);
    }

    void write(const string& result){
        std::lock_guard<std::mutex> lg(mutex);
        outstream << result;
    }

    void flush(){
        outstream.flush();
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

void dispatcher_consumer(ParallelBoundedQueue<Read_Batch>& Q, DispatcherConsumerCallback* cb, LL thread_id){
    write_log("Starting thread " + to_string(thread_id));
    while(true){
        Read_Batch batch  = Q.pop();
        if(batch.data.size() == 0){
            Q.push(batch,0);
            break; // No more work available
        }
        
        Batch_reader br(&batch, 0);
        int64_t read_id = batch.read_id_of_first;
        const char* read;
        LL read_len;
        while(true){
            std::tie(read, read_len) = br.get_next_string();
            if(read_len == 0) break; // End of batch
            cb->callback(read, read_len, read_id);
            read_id++;
        }
    }
    write_log("Thread " + to_string(thread_id) + " done");
}

// Will run characters through fix_char, which at the moment of writing this comment
// upper-cases the character and further the result is not A, C, G or T, changes it to A.
void dispatcher_producer(ParallelBoundedQueue<Read_Batch>& Q, string query_file, int64_t buffer_size){
    // Push work in batches of approximately buffer_size base pairs

    Sequence_Reader fr(query_file, FASTA_MODE);

    Read_Batch rb;
    rb.read_id_of_first = 0;
    int64_t read_id_of_next = 0;

    while(!fr.done()){

        Read_stream input = fr.get_next_query_stream();
        char c;
        bool start = true;
        while(input.getchar(c)){
            rb.data += c;
            rb.read_starts.push_back(start);
            start = false;
        }
        read_id_of_next++;
        
        if(rb.data.size() > buffer_size || fr.done()){
            Q.push(rb,rb.data.size());
            
            // Clear the batch
            rb.data = "";
            rb.read_starts.clear();
            rb.read_id_of_first = read_id_of_next;
        }
    }

    Q.push(rb,0); // Empty batch in the end signifies end of the queue
}

void run_dispatcher(vector<DispatcherConsumerCallback*>& callbacks, string fastafile, LL buffer_size){
    vector<std::thread> threads;
    ParallelBoundedQueue<Read_Batch> Q(buffer_size);

    // Create consumers
    for(int64_t i = 0; i < callbacks.size(); i++){
        threads.push_back(std::thread(dispatcher_consumer,std::ref(Q),callbacks[i], i));
    }

    // Create producer
    threads.push_back(std::thread(dispatcher_producer,std::ref(Q),fastafile,buffer_size));

    for(std::thread& t : threads) t.join();

    for(int64_t i = 0; i < callbacks.size(); i++){
        callbacks[i]->finish();
    }

}
