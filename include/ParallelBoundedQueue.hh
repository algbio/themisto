#pragma once

#include <queue>
#include <thread>
#include <mutex>
#include <iostream>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <utility>
#include <atomic>

using namespace std;

/*

A single producer multiple consumer bounded parallel queue. This implements the solution in Wikipedia: https://en.wikipedia.org/wiki/Monitor_(synchronization)#Condition_variables
(accessed 4.9. 2019). Wikipedia says:

>A classic concurrency problem is that of the bounded producer/consumer, in which there is a queue or ring buffer of tasks with a maximum size, with one or more threads being "producer" threads that add tasks to the queue, and one or more other threads being "consumer" threads that take tasks out of the queue. The queue is assumed to be non–thread-safe itself, and it can be empty, full, or between empty and full. Whenever the queue is full of tasks, then we need the producer threads to block until there is room from consumer threads dequeueing tasks. On the other hand, whenever the queue is empty, then we need the consumer threads to block until more tasks are available due to producer threads adding them.

>As the queue is a concurrent object shared between threads, accesses to it must be made atomic, because the queue can be put into an inconsistent state during the course of the queue access that should never be exposed between threads. Thus, any code that accesses the queue constitutes a critical section that must be synchronized by mutual exclusion. If code and processor instructions in critical sections of code that access the queue could be interleaved by arbitrary context switches between threads on the same processor or by simultaneously-running threads on multiple processors, then there is a risk of exposing inconsistent state and causing race conditions. 

>The solution is to use condition variables. Conceptually a condition variable is a queue of threads, associated with a monitor, on which a thread may wait for some condition to become true. Thus each condition variable c is associated with an assertion P_c. While a thread is waiting on a condition variable, that thread is not considered to occupy the monitor, and so other threads may enter the monitor to change the monitor's state. In most types of monitors, these other threads may signal the condition variable c to indicate that assertion P_c is true in the current state. 

The pseudocode in Wikipedia is:

global volatile RingBuffer queue; // A thread-unsafe ring-buffer of tasks.
global Lock queueLock;  	// A mutex for the ring-buffer of tasks. (Not a spin-lock.)
global CV queueEmptyCV; 	// A condition variable for consumer threads waiting for the queue to 
				// become non-empty.
                        	// Its associated lock is "queueLock".
global CV queueFullCV; 		// A condition variable for producer threads waiting for the queue 
				// to become non-full. Its associated lock is also "queueLock".

// Method representing each producer thread's behavior:
public method producer(){
    while(true){
        task myTask=...; // Producer makes some new task to be added.

        queueLock.acquire(); // Acquire lock for initial predicate check.
        while(queue.isFull()){ // Check if the queue is non-full.
            // Make the threading system atomically release queueLock,
            // enqueue this thread onto queueFullCV, and sleep this thread.
            wait(queueLock, queueFullCV);
            // Then, "wait" automatically re-acquires "queueLock" for re-checking
            // the predicate condition.
        }
        
        // Critical section that requires the queue to be non-full.
        // N.B.: We are holding queueLock.
        queue.enqueue(myTask); // Add the task to the queue.

        // Now the queue is guaranteed to be non-empty, so signal a consumer thread
        // or all consumer threads that might be blocked waiting for the queue to be non-empty:
        signal(queueEmptyCV); -- OR -- notifyAll(queueEmptyCV);
        
        // End of critical sections related to the queue.
        queueLock.release(); // Drop the queue lock until we need it again to add the next task.
    }
}

// Method representing each consumer thread's behavior:
public method consumer(){
    while(true){

        queueLock.acquire(); // Acquire lock for initial predicate check.
        while (queue.isEmpty()){ // Check if the queue is non-empty.
            // Make the threading system atomically release queueLock,
            // enqueue this thread onto queueEmptyCV, and sleep this thread.
            wait(queueLock, queueEmptyCV);
            // Then, "wait" automatically re-acquires "queueLock" for re-checking
            // the predicate condition.
        }
        // Critical section that requires the queue to be non-empty.
        // N.B.: We are holding queueLock.
        myTask=queue.dequeue(); // Take a task off of the queue.
        // Now the queue is guaranteed to be non-full, so signal a producer thread
        // or all producer threads that might be blocked waiting for the queue to be non-full:
        signal(queueFullCV); -- OR -- notifyAll(queueFullCV);

        // End of critical sections related to the queue.
        queueLock.release(); // Drop the queue lock until we need it again to take off the next task.

        doStuff(myTask); // Go off and do something with the task.
    }
}

There is also this codereview.stackoverflow.com post which has this implemented in C++
and nobody has pointed out flaws related to the concurrency:
https://codereview.stackexchange.com/questions/171360/thread-safe-bounded-buffer-fifo-queue-using-c11

*/

template <typename T>
class ParallelBoundedQueue
{

  // Behavior: pop blocks if queue is empty. Push blocks it queue has more load that max_load.

 public:
 
  ParallelBoundedQueue<T>(LL max_load) : current_load(0), max_load(max_load)  {
      assert(max_load > 0);
  }

  T pop(){
    std::unique_lock<std::mutex> lock(queueLock);
    while(queue.empty()) // Check the condition
        queueEmptyCV.wait(lock);

    // Critical section below
    pair<T,LL> item = queue.front(); queue.pop();
    current_load -= item.second;
    queueFullCV.notify_all();
    return item.first;
  } // Lock is released when leaving this function
 
  void push(const T& item, LL load){
    std::unique_lock<std::mutex> lock(queueLock);
    while(current_load > max_load) // Check the condition
        queueFullCV.wait(lock);
    
    // Critical section below
    queue.push({item,load});
    current_load += load;
    queueEmptyCV.notify_all();
  } // Lock is released when leaving the function
 
  private:
  std::queue<pair<T,LL> > queue; // (Element, load) pairs
  std::mutex queueLock;
  std::condition_variable queueEmptyCV;
  std::condition_variable queueFullCV;

  LL current_load;
  const LL max_load;

};
