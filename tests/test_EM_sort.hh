#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <set>
#include <utility>
#include <tuple>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <gtest/gtest.h>
#include "../globals.hh"
#include "../EM_sort.hh"
#include "../test_tools.hh"
#include "setup_tests.hh"

// Returns filename
string generate_line_based_testcase(LL n_lines, LL line_length){
    string outfile = get_temp_file_name("");
    throwing_ofstream out(outfile);
    for(LL i = 0; i < n_lines; i++){
        out << get_random_string(line_length, 2) << "\n";
    }
    out.close();
    return outfile;
}

string sort_lines_stdlib(string infile, const std::function<bool(const char* x, const char* y)>& cmp){
    auto cmp_wrap = [&](const string& x, const string& y){
        return cmp(x.c_str(), y.c_str());
    };
    
    vector<string> lines = get_all_lines(infile);
    sort(lines.begin(), lines.end(), cmp_wrap);
    string outfile = get_temp_file_name("");
    throwing_ofstream out(outfile);
    for(string line : lines) out << line << "\n";
    return outfile;
}

void test_line_sort(string infile, const std::function<bool(const char*, const char*)> & cmp){
    
    string fileA = sort_lines_stdlib(infile, cmp);
    string fileB = temp_file_manager.get_temp_file_name("");
    LL ram = rand() % 10000 + 1;
    LL merge_k = rand() % 6 + 2;
    logger << "Sorting lines with " << ram << " RAM and " << merge_k << "-way merge" << endl;
    EM_sort(infile, fileB, cmp, ram, merge_k, 3, EM_LINES);
    ASSERT_TRUE(files_are_equal(fileA, fileB));

    temp_file_manager.delete_file(fileA);
    temp_file_manager.delete_file(fileB);
    
}

string generate_variable_binary_testcase(LL max_record_len_bytes, LL n_records){
    string outfile = temp_file_manager.get_temp_file_name("");
    throwing_ofstream out(outfile, ios::binary);
    for(LL i = 0; i < n_records; i++){
        LL record_len = max((LL)8, rand() % (max_record_len_bytes + 1));
        //cout << "Add record of length " << record_len << endl;
        write_big_endian_LL(out, record_len);
        for(LL b = 0; b < record_len - 8; b++){
            char byte = static_cast<char>(rand() % 256);
            out.write(&byte, 1);
        }
    }
    return outfile;
}

string generate_constant_binary_testcase(LL record_len, LL n_records){
    string outfile = temp_file_manager.get_temp_file_name("");
    throwing_ofstream out(outfile, ios::binary);
    for(LL i = 0; i < n_records; i++){
        for(LL b = 0; b < record_len; b++){
            char byte = static_cast<char>(rand() % 256);
            out.write(&byte, 1);
        }
    }
    return outfile;
}

string record_to_string(const char* rec){
    stringstream ss;
    LL len = parse_big_endian_LL(rec);
    ss << len << ": ";
    for(LL i = 0; i < len-8; i++){
        int value = static_cast<int>(*reinterpret_cast<const unsigned char*>(rec+8+i));
        ss << value << " ";
    }
    return ss.str();
}

string binary_sort_stdlib(string infile, const std::function<bool(const char*, const char*)> & cmp){
    
    // Using the Block class to help reading in the records
    throwing_ifstream in(infile);
    Variable_binary_block* block = get_next_variable_binary_block(in,(LL)1e16);
    auto cmp_wrap = [&](LL x, LL y){
        return cmp(block->data+x,block->data+y);
    };
    std::sort(block->starts.begin(), block->starts.end(), cmp_wrap);

    // Verify that the sort worked
    for(LL i = 1; i < block->starts.size(); i++){
        // record i must not be strictly less than record i-1
        EXPECT_FALSE(cmp(block->data + block->starts[i], block->data + block->starts[i-1]));
    }

    string outfile = temp_file_manager.get_temp_file_name("");
    throwing_ofstream out(outfile, ios::binary);
    for(LL i = 0; i < block->starts.size(); i++){
        LL length = parse_big_endian_LL(block->data + block->starts[i]);
        out.write(block->data + block->starts[i], length);
    }
    delete block;
    return outfile;
}

string constant_binary_sort_stdlib(string infile, LL rec_len, const std::function<bool(const char*, const char*)> & cmp){
    
    // Using the Block class to help reading in the records
    throwing_ifstream in(infile);
    Constant_binary_block* block = get_next_constant_binary_block(in,(LL)1e16, rec_len);
    auto cmp_wrap = [&](LL x, LL y){
        return cmp(block->data+x,block->data+y);
    };
    std::sort(block->starts.begin(), block->starts.end(), cmp_wrap);

    // Verify that the sort worked
    for(LL i = 1; i < block->starts.size(); i++){
        // record i must not be strictly less than record i-1
        EXPECT_FALSE(cmp(block->data + block->starts[i], block->data + block->starts[i-1]));
    }

    string outfile = temp_file_manager.get_temp_file_name("");
    throwing_ofstream out(outfile, ios::binary);
    for(LL i = 0; i < block->starts.size(); i++){
        out.write(block->data + block->starts[i], rec_len);
    }
    delete block;
    return outfile;
}

void test_variable_binary_sort(string infile, const std::function<bool(const char*, const char*)> & cmp){
    
    string fileA = binary_sort_stdlib(infile, cmp);
    string fileB = temp_file_manager.get_temp_file_name("");
    LL ram = rand() % 10000 + 1;
    LL merge_k = rand() % 6 + 2;
    logger << "Sorting variable binary records with " << ram << " RAM and " << merge_k << "-way merge" << endl;
    EM_sort(infile, fileB, cmp, ram, merge_k, 3, EM_VARIABLE_BINARY);
    ASSERT_TRUE(files_are_equal(fileA, fileB));

    temp_file_manager.delete_file(fileA);
    temp_file_manager.delete_file(fileB);
    
}

void test_constant_binary_sort(string infile, LL record_len, const std::function<bool(const char*, const char*)> & cmp){
    
    string fileA = constant_binary_sort_stdlib(infile, record_len, cmp);
    string fileB = temp_file_manager.get_temp_file_name("");
    LL ram = rand() % 10000 + 1;
    LL merge_k = rand() % 6 + 2;
    logger << "Sorting constant-binary records with " << ram << " RAM and " << merge_k << "-way merge" << endl;
    EM_sort_constant_binary(infile, fileB, cmp, ram, merge_k, record_len, 3);
    ASSERT_TRUE(files_are_equal(fileA, fileB));

    temp_file_manager.delete_file(fileA);
    temp_file_manager.delete_file(fileB);
    
}

TEST(TEST_EM_SORT, variable_binary_sort){
    for(LL max_record_len_bytes = 8; max_record_len_bytes <= 1e6; max_record_len_bytes *= 2){
        // Max record length must be at least 8 because of the long long at the start that tells the length
        for(LL n_records = 1; n_records <= 1e6 && max_record_len_bytes*n_records <= 1e6; n_records *= 2){
            logger << max_record_len_bytes << " " << n_records << endl;
            string infile = generate_variable_binary_testcase(max_record_len_bytes, n_records);
            test_variable_binary_sort(infile, memcmp_variable_binary_records);
        }
    }
}

TEST(TEST_EM_SORT, constant_binary_sort){
    
    for(LL record_len = 1; record_len <= 1e6; record_len *= 2){
        for(LL n_records = 1; n_records <= 1e6 && record_len*n_records <= 1e6; n_records *= 2){
            logger << record_len << " " << n_records << endl;

            auto cmp = [&](const char* x, const char* y){
                return memcmp(x,y,record_len) < 0;
            };
            
            string infile = generate_constant_binary_testcase(record_len, n_records);
            test_constant_binary_sort(infile, record_len, cmp);
        }
    }
}

TEST(TEST_EM_SORT, line_sort){
    for(LL n_lines = 1; n_lines <= 1e6; n_lines *= 2){
        for(LL line_length = 1; line_length <= 1e6 && n_lines*line_length <= 1e6; line_length *= 2){
            logger << n_lines << " " << line_length << endl;
            string infile = generate_line_based_testcase(n_lines, line_length);
            test_line_sort(infile, lex_compare_cstrings);
            test_line_sort(infile, colex_compare_cstrings);
        }
    }

}
