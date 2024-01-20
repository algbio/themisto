#pragma once
#include <cstdint>
#include <vector>
#include <algorithm>

// This is a class where the user can add runs of sets, and in the end
// ask for statistics on the runs of characters in the sets.
class RunTracker{

    private:

    int64_t n_colors;
    std::vector<int64_t> cur_runs; // Index i contains length of current run of color i
    std::vector<int64_t> max_runs; // Index i contain max run length of color i
    std::vector<int64_t> run_counts; // Index i contains number of runs of color i
    std::vector<int64_t> active_colors_list; // Colors that are currently on a run. This list needs to be kept sorted.
    std::vector<bool> seen_colors_bitmap; // A bitmap of all colors that have been seen so far
    std::vector<int64_t> seen_colors_list; // A list of all colors that have been seen so far, i.e. the positions of 1-bits in the seen colors bitmap

    public:

    RunTracker(int64_t n_colors);
    void reset();
    void add_run(const std::vector<int64_t>& colors, int64_t run_length);
    const std::vector<int64_t>& get_run_counts() const;
    const std::vector<int64_t>& get_max_runs() const;

};