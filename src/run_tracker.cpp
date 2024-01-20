#include <cstdint>
#include <vector>
#include <algorithm>
#include <iostream>
#include "run_tracker.hh"

RunTracker::RunTracker(int64_t n_colors) : n_colors(n_colors) {
    cur_runs.resize(n_colors);
    max_runs.resize(n_colors);
    run_counts.resize(n_colors);
    seen_colors_bitmap.resize(n_colors);
}

void RunTracker::reset(){
    for(int64_t x : seen_colors_list){
        cur_runs[x] = 0;
        max_runs[x] = 0;
        run_counts[x] = 0;
        seen_colors_bitmap[x] = false;
    }
    seen_colors_list.clear();
    active_colors_list.clear();
}

void RunTracker::add_run(const std::vector<int64_t>& colors, int64_t run_length){
    std::vector<int64_t> sorted_colors = colors; // TODO: can we assume that the colors are already sorted?
    std::sort(sorted_colors.begin(), sorted_colors.end());

    // Make note of any colors that we have not yet seen
    for(int64_t x : colors){
        if(!seen_colors_bitmap[x]){
            seen_colors_bitmap[x] = true;
            seen_colors_list.push_back(x);
        }
    }

    // Scan sorted lists A (=active_colors_list) and B (=colors) and classify elements into three categories:
    // 1. Elements in A but not in B
    // 2. Elements in B but not in A
    // 3. Elements in both A and B
    int64_t i = 0, j = 0;
    while(i < active_colors_list.size() || j < colors.size()){
        if(i != active_colors_list.size() && (j == colors.size() || active_colors_list[i] < colors[j])){
            // Was in previous run but not in this one
            cur_runs[active_colors_list[i]] = 0; // Run end
            i++; 
        }
        else if(j != colors.size() && (i == active_colors_list.size() || active_colors_list[i] > colors[j])){
            // Was not in previous run but is in this one
            cur_runs[colors[j]] = run_length; // Run start
            max_runs[colors[j]] = std::max(max_runs[colors[j]], run_length);
            run_counts[colors[j]]++;
            j++;
        }
        else{
            // Is in both
            cur_runs[colors[j]] += run_length; // Run continues
            max_runs[colors[j]] = std::max(max_runs[colors[j]], cur_runs[colors[j]]);
            i++; j++;
        }
    }
    active_colors_list = colors;
}

const std::vector<int64_t>& RunTracker::get_run_counts() const{
    return run_counts;
}

const std::vector<int64_t>& RunTracker::get_max_runs() const{
    return max_runs;
}

const std::vector<int64_t>& RunTracker::get_seen_colors() const{
    return seen_colors_list;
}

int test_run_tracker(){

    RunTracker rt(4);
    rt.add_run({0,1,2}, 12);
    rt.add_run({0,1,  3}, 14);
    rt.add_run({0,  2,3}, 15);
    rt.reset();

    rt.add_run({0,1,2}, 3);
    rt.add_run({0,1,  3}, 4);
    rt.add_run({0,  2,3}, 5);
    rt.add_run({},1);
    rt.add_run({0,1,2,3}, 6);
    rt.add_run({  1,2 }, 6);
    rt.add_run({    2,3}, 5);
    rt.add_run({0    ,3}, 4);
    rt.add_run({0,1   }, 3);
    rt.add_run({  1,2 }, 2);

    std::cout << "--" << std::endl;
    std::vector<int64_t> max_runs = rt.get_max_runs();
    std::vector<int64_t> run_counts = rt.get_run_counts();
    std::vector<int64_t> true_max_runs = {12,12,17,9};
    std::vector<int64_t> true_run_counts = {3,3,4,3};

    for(auto x : run_counts) std::cout << x << " "; 
    std::cout << std::endl;
    for(auto x : true_run_counts) std::cout << x << " "; 
    std::cout << std::endl;

    for(auto x : max_runs) std::cout << x << " ";
    std::cout << std::endl;
    for(auto x : true_max_runs) std::cout << x << " ";
    std::cout << std::endl;
}