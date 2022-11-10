#include <string>
#include <vector>
#include "sbwt/SeqIO.hh"
#include "globals.hh"
#include "zstr.hpp"
#include "WorkDispatcher.hh"
#include "stdlib_printing.hh"

using namespace std;
using namespace sbwt;

void dispatcher_consumer(ParallelBoundedQueue<ReadBatch*>& Q, DispatcherConsumerCallback* cb, int64_t thread_id){
    write_log("Starting thread " + to_string(thread_id), LogLevel::MINOR);
    while(true){
        ReadBatch* batch  = Q.pop();
        if(batch->data.size() == 0){
            Q.push(batch,0);
            break; // No more work available
        }
        ReadBatchIterator rbi(batch, 0);
        int64_t read_id = batch->firstReadID;
        void* metadata = nullptr;
        const char* read;
        int64_t read_len;
        while(true){
            std::tie(read, read_len, metadata) = rbi.getNextRead();
            if(read_len == 0) break; // End of batch
            cb->callback(read, read_len, read_id, metadata);
            read_id++;
        }
        delete batch; // Allocated by producer
    }
    write_log("Thread " + to_string(thread_id) + " done", LogLevel::MINOR);
}
