#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include <cassert>
#include <set>
#include <unordered_map>
#include <map>
#include <gtest/gtest.h>
#include <cassert>
#include "test_tools.hh"
#include "Succinct_Prefix_Sum.hh"

using namespace sbwt;

vector<int64_t> get_sparse_colorset(){
    vector<int64_t> v = {4, 1534, 4003, 8903};
    return v;
}

vector<int64_t> get_dense_colorset(int64_t gap, int64_t total_length){
    vector<int64_t> v;
    for(int64_t i = 0; i < total_length; i += gap){
        v.push_back(i);
    }
    return v;
}

TEST(NEW_NEW_COLORING_TEST, storage){

    New_Color_Set_Storage css;
    vector<vector<int64_t> > sets = {get_sparse_colorset(), 
                                     get_dense_colorset(1,1000), 
                                     get_sparse_colorset(), 
                                     get_sparse_colorset(), 
                                     get_dense_colorset(2,1000),
                                     {1,5,7,8},
                                     {0,1,2,3,4,5,6},
                                     get_dense_colorset(3,1000)};

    for(const vector<int64_t>& set : sets)
        css.add_set(set);
    css.prepare_for_queries();

    vector<Color_Set_View> retrieved_views = css.get_all_sets();
    for(int64_t i = 0; i < retrieved_views.size(); i++){
        for(int64_t x : retrieved_views[i].get_colors_as_vector()) cout << x << " "; 
        cout << endl;

        ASSERT_EQ(retrieved_views[i].get_colors_as_vector(), sets[i]);
    }
    
}

TEST(NEW_NEW_COLORING_TEST, prefix_sums){
    vector<int64_t> v = {0,2,5,10,0,0,0,2,0,4};
    Succinct_Prefix_Sums sps;
    for(int64_t x : v) sps.add(x);
    sps.finish_building();

    vector<int64_t> sums = {0};
    for(int64_t x : v) sums.push_back(sums.back() + x);

    vector<int64_t> test;
    for(int64_t i = 0; i <= v.size(); i++)
        test.push_back(sps.sum(i));

    ASSERT_EQ(sums, test);

    Succinct_Prefix_Sums sps2 = to_disk_and_back(sps);

    vector<int64_t> test2;
    for(int64_t i = 0; i <= v.size(); i++)
        test2.push_back(sps2.sum(i));

    ASSERT_EQ(sums, test2);

}