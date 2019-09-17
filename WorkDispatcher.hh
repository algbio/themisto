#pragma once
#include <string>
#include <vector>
#include "ParallelQueue.hh"
#include "input_reading.hh"
#include "globals.hh"

using namespace std;

class ParallelOutputWriter{
    public:

    string outfile;
    ofstream outstream;
    std::mutex mutex;

    ParallelOutputWriter(string outfile) : outfile(outfile){
        check_writable(outfile);
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

void dispatcher_consumer(ParallelQueue<Read_Batch>& Q, DispatcherConsumerCallback* cb){
    while(true){
        Read_Batch batch  = Q.pop();
        if(batch.data.size() == 0){
            Q.push(batch);
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
}

void dispatcher_producer(ParallelQueue<Read_Batch>& Q, string query_file, int64_t buffer_size){
    // Todo: Q should have limited size?
    // Push work in batches of approximately buffer_size base pairs

    FASTA_reader fr(query_file);

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
            Q.push(rb);
            
            // Clear the batch
            rb.data = "";
            rb.read_starts.clear();
            rb.read_id_of_first = read_id_of_next;
        }
    }

    Q.push(rb); // Empty batch in the end signifies end of the queue
}

void run_dispatcher(vector<DispatcherConsumerCallback*>& callbacks, string fastafile, LL buffer_size){
    vector<std::thread> threads;
    ParallelQueue<Read_Batch> Q;

    // Create consumers
    for(int64_t i = 0; i < callbacks.size(); i++){
        threads.push_back(std::thread(dispatcher_consumer,std::ref(Q),callbacks[i]));
    }

    // Create producer
    threads.push_back(std::thread(dispatcher_producer,std::ref(Q),fastafile,buffer_size));

    for(std::thread& t : threads) t.join();

    for(int64_t i = 0; i < callbacks.size(); i++){
        callbacks[i]->finish();
    }

}
