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

    Color_Set_Storage<Color_Set> css;
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

    // Check that we can get back the same color sets as what we put in
    vector<Color_Set::view_t> retrieved_views = css.get_all_sets();
    for(int64_t i = 0; i < retrieved_views.size(); i++){
        ASSERT_EQ(retrieved_views[i].get_colors_as_vector(), sets[i]);
        ASSERT_FALSE(retrieved_views[i].empty());
        ASSERT_EQ(retrieved_views[i].size(), sets[i].size());

        // Check contains()

        vector<bool> contains_ref(10000);
        for(int64_t x : sets[i]) contains_ref[x] = true;    

        vector<bool> contains_check(10000);
        for(int64_t j = 0; j < contains_check.size(); j++) contains_check[j] = retrieved_views[i].contains(j);

        ASSERT_EQ(contains_check, contains_ref);

        // Test constructing a color set object out of a view
        Color_Set cs(retrieved_views[i]);
        ASSERT_EQ(cs.get_colors_as_vector(), retrieved_views[i].get_colors_as_vector());
        ASSERT_EQ(cs.empty(), retrieved_views[i].empty());
        ASSERT_EQ(cs.size(), retrieved_views[i].size());
        for(int64_t j = 0; j < 10000; j++){
            ASSERT_EQ(cs.contains(j), retrieved_views[i].contains(j));            
        }

        // Test copy
        Color_Set cs2(cs);
        ASSERT_EQ(cs.get_colors_as_vector(), cs2.get_colors_as_vector());
        ASSERT_EQ(cs.empty(), cs2.empty());
        ASSERT_EQ(cs.size(), cs2.size());
        for(int64_t j = 0; j < 10000; j++){
            ASSERT_EQ(cs.contains(j), cs2.contains(j));            
        }
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