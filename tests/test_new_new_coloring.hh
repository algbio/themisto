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