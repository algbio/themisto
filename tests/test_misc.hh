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
#include "../stdlib_printing.hh"
#include "../globals.hh"
#include "setup_tests.hh"
#include "test_tools.hh"
#include "Themisto.hh"
#include "Argv.hh"
#include "commands.hh"
#include <cassert>

TEST(MISC_TEST, delete_non_ACGT){

    get_temp_file_manager().set_dir("temp");

    string fasta_data = ">\nNNAGTATTGNNCTGTAGXGTCAGTGCTACGTACTNN\n>\nAGTTCTNNTAGTGCNNTNXNAGCCA\n";
    //                      ^^       ^^      ^                ^^           ^^      ^^ ^^^
    string color_data = "0\n1\n";

    vector<string> correct_out_seqs = {"AGTATTG", "CTGTAG", "GTCAGTGCTACGTACT", "AGTTCT", "TAGTGC", "T", "AGCCA"};
    vector<LL> correct_out_colors = {0,0,0,1,1,1,1};

    // Write to files

    string fasta1 = get_temp_file_manager().create_filename("fasta1");
    string colors1 = get_temp_file_manager().create_filename("colors1");
    
    throwing_ofstream fasta1_outstream(fasta1);
    fasta1_outstream << fasta_data;
    fasta1_outstream.close();

    throwing_ofstream colors1_outstream(colors1);
    colors1_outstream << color_data;
    colors1_outstream.close();

    // Split

    string fasta2, colors2;
    std::tie(fasta2, colors2) = split_all_seqs_at_non_ACGT(fasta1, "fasta", colors1);

    // Verify
    vector<LL> out_colors;
    string line;
    ifstream colors_stream(colors2);
    while(getline(colors_stream, line))
        out_colors.push_back(stoll(line));
    vector<string> out_seqs;
    Sequence_Reader sr(fasta2, FASTA_MODE);
    while(!sr.done()) out_seqs.push_back(sr.get_next_query_stream().get_all());

    //log << correct_out_seqs << endl << out_seqs << endl << correct_out_colors << endl << out_colors << endl;
    logger << "test" << 2 << out_seqs << std::endl<char, std::char_traits<char>>;
    ASSERT_EQ(correct_out_seqs, out_seqs);
    ASSERT_EQ(correct_out_colors, out_colors);

}

TEST(MISC_TEST, string_to_integer_safe){
    try{
        string_to_integer_safe("1 2 3");
        FAIL(); // Should have thrown an exception
    } catch(...){}
    try{
        string_to_integer_safe("  12 33 ");
        FAIL(); // Should have thrown an exception
    } catch(...){}
    try{
        string_to_integer_safe("  12X33 ");
        FAIL(); // Should have thrown an exception
    } catch(...){}
    try{
        string_to_integer_safe("  123XX ");
        FAIL(); // Should have thrown an exception
    } catch(...){}
    try{
        string_to_integer_safe("   ");
        FAIL(); // Should have thrown an exception
    } catch(...){}

    ASSERT_EQ(string_to_integer_safe("  \n\t  \r 1234567890\n  \r\n"), 1234567890);
}
