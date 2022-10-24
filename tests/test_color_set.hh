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
#include "hybrid_color_set.hh"
#include "Fixed_Width_Int_Color_Set.hh"
#include "bit_magic_color_set.hh"
#include "new_new_coloring.hh"

using namespace sbwt;

vector<int64_t> get_sparse_example(){
    vector<int64_t> v = {4, 1534, 4003, 8903};
    return v;
}

vector<int64_t> get_dense_example(int64_t gap, int64_t total_length){
    vector<int64_t> v;
    for(int64_t i = 0; i < total_length; i += gap){
        v.push_back(i);
    }
    return v;
}

template<typename color_set_t> requires Color_Set_Interface<color_set_t>
void test_sparse_color_set(){
    vector<int64_t> v = get_sparse_example();
    color_set_t cs(v);

    vector<int64_t> v2 = cs.get_colors_as_vector();

    ASSERT_EQ(v,v2);

    int64_t n = *std::max_element(v.begin(), v.end());
    for(int64_t x = 0; x <= n + 10; x++){ // +10: test over the boundary too.
        bool found = std::find(v.begin(), v.end(), x) != v.end();
        ASSERT_EQ(cs.contains(x), found);
    }
}

TEST(TEST_COLOR_SET, sparse){
    test_sparse_color_set<Roaring_Color_Set>();
    test_sparse_color_set<Bit_Magic_Color_Set>();
    test_sparse_color_set<Color_Set>();
}

template<typename color_set_t> requires Color_Set_Interface<color_set_t>
void test_dense_color_set(){
    vector<int64_t> v = get_dense_example(3, 1000);
    color_set_t cs(v);

    vector<int64_t> v2 = cs.get_colors_as_vector();

    ASSERT_EQ(v,v2);

    for(int64_t i = 0; i < 1000 + 10; i++){
        ASSERT_EQ(cs.contains(i), (i < 1000 && i % 3 == 0));
    }
}

TEST(TEST_COLOR_SET, dense){
    test_dense_color_set<Roaring_Color_Set>();
    test_dense_color_set<Bit_Magic_Color_Set>();
    test_dense_color_set<Color_Set>();
}

template<typename color_set_t> requires Color_Set_Interface<color_set_t>
void test_sparse_vs_sparse(){
    vector<int64_t> v1 = {4, 1534, 4003, 8903};
    vector<int64_t> v2 = {4, 2000, 4003, 5000};

    color_set_t c1(v1);
    color_set_t c2(v2);
    color_set_t c12(c1);
    cout << "==" << endl;
    print(c1.get_colors_as_vector());
    print(c2.get_colors_as_vector());
    print(c12.get_colors_as_vector());
    cout << "==" << endl;
    c12.intersection(c2);

    vector<int64_t> v12_inter = c12.get_colors_as_vector();
    vector<int64_t> correct_inter = {4,4003};
    ASSERT_EQ(v12_inter, correct_inter);

    color_set_t c12_union(c1);
    c12_union.do_union(c2);
    vector<int64_t> v12_union = c12_union.get_colors_as_vector();
    vector<int64_t> correct_union = {4, 1534, 2000, 4003, 5000, 8903};
    ASSERT_EQ(v12_union, correct_union);
}
TEST(TEST_COLOR_SET, sparse_vs_sparse){
    test_sparse_vs_sparse<Roaring_Color_Set>();
    test_sparse_vs_sparse<Bit_Magic_Color_Set>();
    test_sparse_vs_sparse<Color_Set>();
}

template<typename color_set_t> requires Color_Set_Interface<color_set_t>
void test_dense_vs_dense(){
    vector<int64_t> v1 = get_dense_example(2, 1000); // Multiples of 2
    vector<int64_t> v2 = get_dense_example(3, 1000); // Multiples of 3

    color_set_t c1(v1);
    color_set_t c2(v2);
    color_set_t c12(c1);
    c12.intersection(c2);

    vector<int64_t> v12_inter = c12.get_colors_as_vector();
    vector<int64_t> correct_inter = get_dense_example(6, 1000); // 6 = lcm(2,3)
    ASSERT_EQ(v12_inter, correct_inter);

    color_set_t c12_union(c1);
    c12_union.do_union(c2);
    vector<int64_t> v12_union = c12_union.get_colors_as_vector();
    vector<int64_t> correct_union;
    for(int64_t i = 0; i < 1000; i++){
        if(i % 2 == 0 || i % 3 == 0) correct_union.push_back(i);
    }
    ASSERT_EQ(v12_union, correct_union);
}

TEST(TEST_COLOR_SET, dense_vs_dense){
    test_dense_vs_dense<Roaring_Color_Set>();
    test_dense_vs_dense<Bit_Magic_Color_Set>();
    test_dense_vs_dense<Color_Set>();
}

template<typename color_set_t> requires Color_Set_Interface<color_set_t>
void test_dense_vs_sparse(){
    vector<int64_t> v1 = get_dense_example(3, 10000); // Multiples of 3
    vector<int64_t> v2 = {3, 4, 5, 3000, 6001, 9999};

    color_set_t c1(v1);
    color_set_t c2(v2);
    color_set_t c12(c1);
    c12.intersection(c2);

    vector<int64_t> v12_inter = c12.get_colors_as_vector();
    vector<int64_t> correct_inter = {3, 3000, 9999};
    ASSERT_EQ(v12_inter, correct_inter);

    color_set_t c12_union(c1);
    c12_union.do_union(c2);
    vector<int64_t> v12_union = c12_union.get_colors_as_vector();
    vector<int64_t> correct_union;
    for(int64_t i = 0; i < 10000; i++){
        if(i % 3 == 0 || std::find(v2.begin(), v2.end(), i) != v2.end()) correct_union.push_back(i);
    }

    ASSERT_EQ(v12_union, correct_union);
}

TEST(TEST_COLOR_SET, dense_vs_sparse){
    test_dense_vs_sparse<Roaring_Color_Set>();
    test_dense_vs_sparse<Bit_Magic_Color_Set>();
    test_dense_vs_sparse<Color_Set>();
}


template<typename color_set_t> requires Color_Set_Interface<color_set_t>
void test_sparse_vs_dense(){
    vector<int64_t> v1 = {3, 4, 5, 3000, 6001, 9999};
    vector<int64_t> v2 = get_dense_example(3, 10000); // Multiples of 3

    color_set_t c1(v1);
    color_set_t c2(v2);
    color_set_t c12(c1);
    c12.intersection(c2);

    vector<int64_t> v12_inter = c12.get_colors_as_vector();
    vector<int64_t> correct_inter = {3, 3000, 9999};
    ASSERT_EQ(v12_inter, correct_inter);

    color_set_t c12_union(c1);
    c12_union.do_union(c2);
    vector<int64_t> v12_union = c12_union.get_colors_as_vector();
    vector<int64_t> correct_union;
    for(int64_t i = 0; i < 10000; i++){
        if(i % 3 == 0 || std::find(v2.begin(), v2.end(), i) != v2.end()) correct_union.push_back(i);
    }

    ASSERT_EQ(v12_union, correct_union);
}

TEST(TEST_COLOR_SET, sparse_vs_dense){
    test_sparse_vs_dense<Roaring_Color_Set>();
    test_sparse_vs_dense<Bit_Magic_Color_Set>();
    test_sparse_vs_dense<Color_Set>();
}

template<typename color_set_t> requires Color_Set_Interface<color_set_t>
void test_empty_color_set(){
    vector<int64_t> v = get_sparse_example();

    color_set_t c1(v);
    color_set_t c2;

    ASSERT_FALSE(c1.empty());
    ASSERT_TRUE(c2.empty());
}

TEST(TEST_COLOR_SET, empty){
    test_empty_color_set<Roaring_Color_Set>();
    test_empty_color_set<Bit_Magic_Color_Set>();
    test_empty_color_set<Color_Set>();
}

template<typename color_set_t> requires Color_Set_Interface<color_set_t>
void test_dense_color_set_serialization(){
    vector<int64_t> v = get_dense_example(3, 10000); // Multiples of 3
    color_set_t c(v);
    color_set_t c2 = to_disk_and_back(c);
    ASSERT_EQ(c2.get_colors_as_vector(), v);
}

TEST(TEST_COLOR_SET, dense_serialization){
    test_dense_color_set_serialization<Roaring_Color_Set>();
    test_dense_color_set_serialization<Bit_Magic_Color_Set>();
    // The Newnew color set does not have serialization for individual color sets
}

template<typename color_set_t> requires Color_Set_Interface<color_set_t>
void test_sparse_color_set_serialization(){
    vector<int64_t> v = get_sparse_example();
    color_set_t c(v);
    color_set_t c2 = to_disk_and_back(c);
    ASSERT_EQ(c2.get_colors_as_vector(), v);
}

TEST(TEST_COLOR_SET, sparse_serialization){
    test_sparse_color_set_serialization<Roaring_Color_Set>();
    test_sparse_color_set_serialization<Bit_Magic_Color_Set>();
    // The Newnew color set does not have serialization for individual color sets
}
