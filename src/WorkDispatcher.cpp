#pragma once
#include <string>
#include <vector>
#include "ParallelBoundedQueue.hh"
#include "input_reading.hh"
#include "globals.hh"
#include "zstr.hpp"
#include "WorkDispatcher.hh"

using namespace std;

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
void dispatcher_producer(ParallelBoundedQueue<Read_Batch>& Q, Sequence_Reader& sr, int64_t buffer_size){
    // Push work in batches of approximately buffer_size base pairs

    Read_Batch rb;
    rb.read_id_of_first = 0;
    int64_t read_id_of_next = 0;

    while(!sr.done()){

        Read_stream input = sr.get_next_query_stream();
        char c;
        bool start = true;
        while(input.getchar(c)){
            rb.data += c;
            rb.read_starts.push_back(start);
            start = false;
        }
        read_id_of_next++;
        
        if(rb.data.size() > buffer_size || sr.done()){
            Q.push(rb,rb.data.size());
            
            // Clear the batch
            rb.data = "";
            rb.read_starts.clear();
            rb.read_id_of_first = read_id_of_next;
        }
    }

    Q.push(rb,0); // Empty batch in the end signifies end of the queue
}

void run_dispatcher(vector<DispatcherConsumerCallback*>& callbacks, Sequence_Reader& sr, LL buffer_size){
    vector<std::thread> threads;
    ParallelBoundedQueue<Read_Batch> Q(buffer_size);

    // Create consumers
    for(int64_t i = 0; i < callbacks.size(); i++){
        threads.push_back(std::thread(dispatcher_consumer,std::ref(Q),callbacks[i], i));
    }

    // Create producer
    threads.push_back(std::thread(dispatcher_producer,std::ref(Q),std::ref(sr),buffer_size));

    for(std::thread& t : threads) t.join();

    for(int64_t i = 0; i < callbacks.size(); i++){
        callbacks[i]->finish();
    }

}
