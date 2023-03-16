#include <vector>
#include "ParallelBoundedQueue.hh"
#include <mutex>
#include <cmath>
#include <optional>
#include <numeric>

using namespace std;

class ExampleWorkerThread{

    private:

    vector<int64_t> results;
    std::mutex* mutex;
    int64_t n_calls_between_flushes;

    public:

    struct WorkItem{
        vector<int64_t> data;
    };

    // Takes in the mutex that protects flush_results
    ExampleWorkerThread(std::mutex* mutex, int64_t n_calls_between_flushes) : n_calls_between_flushes(n_calls_between_flushes){
        this->mutex = mutex;
    }

    // This function should only use local variables and no shared state
    void run(sbwt::ParallelBoundedQueue<std::optional<WorkItem>>& Q){
        int64_t calls_to_next_flush = n_calls_between_flushes;
        while(std::optional<WorkItem> item = Q.pop()){
            
            // Process the work item
            int64_t sum = 0;
            for(int64_t i = 0; i < 1000000; i++){
                sum += std::accumulate(item->data.begin(), item->data.end(), 0);
            }
            results.push_back(sum);

            // Flush results if needed
            calls_to_next_flush--;
            if(calls_to_next_flush == 0){
                calls_to_next_flush = n_calls_between_flushes;
                flush_results();
            }
        }

        flush_results();
    }

    // This function may use shared state.
    void flush_results(){
        std::lock_guard<std::mutex> lock(*mutex);

        for(int64_t x : results) cout << x << "\n";
        
        results.clear();
    }

};

template<typename worker_t>
class ThreadPool{

    typedef typename worker_t::WorkItem work_item_t;

    int64_t next_worker_id;
    vector<worker_t> workers;
    vector<std::thread> threads;
    sbwt::ParallelBoundedQueue<std::optional<work_item_t>> work_queue;
    std::mutex mutex;

    public:

    ThreadPool(int64_t n_workers, int64_t max_work_queue_load, int64_t n_calls_between_output_flushes) : next_worker_id(0), work_queue(max_work_queue_load){
        for(int64_t i = 0; i < n_workers; i++){
            workers.push_back(worker_t(&mutex, n_calls_between_output_flushes));
            threads.push_back(
                std::thread([this, i]{
                    this->workers[i].run(this->work_queue);
                })
            );
        }
    }

    void add_work(work_item_t input, int64_t load){
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