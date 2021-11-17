#include <string>
#include <vector>
#include "ParallelBoundedQueue.hh"
#include "input_reading.hh"
#include "globals.hh"
#include "zstr.hpp"
#include "ReadBatch.hh"
#include "WorkDispatcher.hh"

using namespace std;

void dispatcher_consumer(ParallelBoundedQueue<ReadBatch*>& Q, DispatcherConsumerCallback* cb, LL thread_id){
    write_log("Starting thread " + to_string(thread_id));
    while(true){
        ReadBatch* batch  = Q.pop();
        if(batch->data.size() == 0){
            Q.push(batch,0);
            break; // No more work available
        }
        
        ReadBatchIterator rbi(batch, 0);
        int64_t read_id = batch->firstReadID;
        const char* read;
        LL read_len;
        while(true){
            std::tie(read, read_len) = rbi.getNextRead();
            if(read_len == 0) break; // End of batch
            cb->callback(read, read_len, read_id);
            read_id++;
        }
        delete batch; // Allocated by producer
    }
    write_log("Thread " + to_string(thread_id) + " done");
}

void dispatcher_producer(ParallelBoundedQueue<ReadBatch*>& Q, Sequence_Reader_Buffered& sr, int64_t batch_size){
    // Push work in batches of approximately buffer_size base pairs

    LL read_id = 0;
    ReadBatch* batch = new ReadBatch(); // Deleted by consumer

    while(true){
        // Create a new read batch and push it to the work queue

        LL len = sr.get_next_read_to_buffer();
        if(len == 0){
            // All reads read. Push the current batch if it's non-empty and quit
            if(batch->data.size() > 0) {
                Q.push(batch, batch->data.size());
                batch = new ReadBatch(); // Clear
            }
            break; // quit
        } else{
            // Add read to batch
            if(batch->data.size() == 0) batch->firstReadID = read_id++;
            batch->readStarts.push_back(batch->data.size());   
            for(LL i = 0; i < len; i++)
                batch->data.push_back(sr.read_buf[i]);
            if(batch->data.size() >= batch_size){
                Q.push(batch, batch->data.size());
                batch = new ReadBatch(); // Clear
            }
        }
    }
    
    Q.push(batch,0); // Empty batch in the end signifies end of the queue
}

void run_dispatcher(vector<DispatcherConsumerCallback*>& callbacks, Sequence_Reader_Buffered& sr, LL buffer_size){
    vector<std::thread> threads;
    ParallelBoundedQueue<ReadBatch*> Q(buffer_size);

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
