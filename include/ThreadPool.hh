#pragma once

#include <vector>
#include "sbwt/EM_sort/ParallelBoundedQueue.hh"
#include <mutex>
#include <cmath>
#include <optional>
#include <numeric>

using namespace std;

// This is like the parallel bounded queue in SBWT, but uses move-semantics
// inside the queue to avoid copies.
template <typename T>
class ThreadPoolParallelBoundedQueue
{

  // Behavior: pop blocks if queue is empty. Push blocks it queue has more load that max_load.

 public:
 
  ThreadPoolParallelBoundedQueue(int64_t max_load) : current_load(0), max_load(max_load)  {
      assert(max_load > 0);
  }

  T pop(){
    std::unique_lock<std::mutex> lock(queueLock);
    while(queue.empty()) // Check the condition
        queueEmptyCV.wait(lock);

    // Critical section below
    pair<T,int64_t> item = std::move(queue.front()); queue.pop();
    current_load -= item.second;
    queueFullCV.notify_all();
    return std::move(item.first);
  } // Lock is released when leaving this function
 
  void push(T item, int64_t load){
    std::unique_lock<std::mutex> lock(queueLock);
    while(current_load > max_load) // Check the condition
        queueFullCV.wait(lock);
    
    // Critical section below
    queue.push(move(make_pair(move(item),load)));
    current_load += load;
    queueEmptyCV.notify_all();
  } // Lock is released when leaving the function
 
  private:
  std::queue<pair<T,int64_t> > queue; // (Element, load) pairs
  std::mutex queueLock;
  std::condition_variable queueEmptyCV;
  std::condition_variable queueFullCV;

  int64_t current_load;
  const int64_t max_load;

};

// Worker threads should inherit from this class.
// For the ThreadPool to work, the derived worker types must
// implement the virtual methods `process_work_item` and `critical_section`.
template<typename work_item_t>
class BaseWorkerThread{

    private:

    std::mutex* critical_section_mutex;

    public:

    // Called by the ThreadPool
    void set_critical_section_mutex(std::mutex* critical_section_mutex){
        this->critical_section_mutex = critical_section_mutex;
    }

    // This function should only use local variables and no unprotected shared state
    virtual void process_work_item(work_item_t item) = 0;

    // This function is called after every work item. It is called so
    // that only one thread at a time is executing the function.
    virtual void critical_section() = 0;

    
    void run(ThreadPoolParallelBoundedQueue<std::optional<work_item_t>>& Q){
        while(std::optional<work_item_t> item = move(Q.pop())){
            // Process the work item
            process_work_item(std::move(*item));
            std::lock_guard<std::mutex> lock(*critical_section_mutex);
            critical_section();
        }
    }

    virtual ~BaseWorkerThread() = default;

};

template<typename worker_t, typename work_item_t>
class ThreadPool{

    vector<std::thread> threads;
    ThreadPoolParallelBoundedQueue<std::optional<work_item_t>> work_queue;
    std::mutex critical_section_mutex;

    public:

    ThreadPool(vector<worker_t*>& workers, int64_t max_work_queue_load) : work_queue(max_work_queue_load){
        for(worker_t* worker : workers){
            worker->set_critical_section_mutex(&critical_section_mutex);
            threads.push_back(
                std::thread([worker, this]{
                    worker->run(this->work_queue);
                })
            );
        }
    }

    void add_work(work_item_t input, int64_t load){
        work_queue.push(std::move(input), load);
    }

    // Waits until the work queue is empty and all threads are finished
    void join_threads(){
        // Add null work items to signify the end of the queue
        for(int64_t i = 0; i < threads.size(); i++)
            work_queue.push(std::nullopt, 0);
        for(auto& t : threads) t.join();
    }

};