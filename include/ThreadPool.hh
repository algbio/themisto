#pragma once

#include <vector>
#include "sbwt/EM_sort/ParallelBoundedQueue.hh"
#include <mutex>
#include <cmath>
#include <optional>
#include <numeric>

using namespace std;
template<typename work_item_t>
class BaseWorkerThread{

    private:

    std::mutex* mutex;
    int64_t n_calls_between_flushes;

    public:


    // Takes in the mutex that protects flush_results
    BaseWorkerThread(std::mutex* mutex, int64_t n_calls_between_flushes) : n_calls_between_flushes(n_calls_between_flushes){
        this->mutex = mutex;
    }

    // This function should only use local variables and no shared state
    virtual void process_work_item(work_item_t item) = 0;

    // This function is protected with a lock and may use shared state.
    virtual void flush_results() = 0;

    
    void run(sbwt::ParallelBoundedQueue<std::optional<work_item_t>>& Q){

        int64_t calls_to_next_flush = n_calls_between_flushes;
        while(std::optional<work_item_t> item = Q.pop()){
            
            // Process the work item
            process_work_item(*item);

            // Flush results if needed
            calls_to_next_flush--;
            if(calls_to_next_flush == 0){
                calls_to_next_flush = n_calls_between_flushes;
                std::lock_guard<std::mutex> lock(*mutex);
                flush_results();
            }
        }

        std::lock_guard<std::mutex> lock(*mutex);
        flush_results();
    }

};


struct ExampleWorkItem{
    vector<int64_t> data;
};

class ExampleWorkerThread : public BaseWorkerThread<ExampleWorkItem>{

    private:

    vector<int64_t> results;

    public:

    typedef ExampleWorkItem work_item_t;

    ExampleWorkerThread(std::mutex* mutex, int64_t n_calls_between_flushes) : BaseWorkerThread(mutex, n_calls_between_flushes) {}

    // This function should only use local variables and no shared state
    virtual void process_work_item(ExampleWorkItem item){
        // Do some computation on the item
        int64_t sum = 0;
        while(item.data.size() > 0){
            sum += item.data[0];
            item.data.erase(item.data.begin());
        }
        results.push_back(sum);
    }

    // This function is protected with a lock and may use shared state.
    virtual void flush_results(){
        for(int64_t x : results) cout << x << "\n";
        results.clear();
    }

};

template<typename worker_t>
class ThreadPool{

    int64_t next_worker_id;
    vector<unique_ptr<BaseWorkerThread<typename worker_t::work_item_t>>> workers;
    vector<std::thread> threads;
    sbwt::ParallelBoundedQueue<std::optional<typename worker_t::work_item_t>> work_queue;
    std::mutex mutex;

    public:

    ThreadPool(int64_t n_workers, int64_t max_work_queue_load, int64_t n_calls_between_output_flushes) : next_worker_id(0), work_queue(max_work_queue_load){
        for(int64_t i = 0; i < n_workers; i++){
            unique_ptr<BaseWorkerThread<typename worker_t::work_item_t>> worker = make_unique<ExampleWorkerThread>(&mutex, n_calls_between_output_flushes);
            workers.push_back(std::move(worker));
            threads.push_back(
                std::thread([this, i]{
                    this->workers[i]->run(this->work_queue);
                })
            );
        }
    }

    void add_work(typename worker_t::work_item_t input, int64_t load){
        work_queue.push(input, load);
    }

    // Waits until the work queue is empty and all threads are finished
    void join_threads(){
        // Add null work items to signify the end of the queue
        for(int64_t i = 0; i < workers.size(); i++)
            work_queue.push(std::nullopt, 0);
        for(auto& t : threads) t.join();
    }

};