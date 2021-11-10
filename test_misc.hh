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
#include "globals.hh"
#include <cassert>

bool misc_test_delete_non_ACGT(){

    string fasta_data = ">\nNNAGTATTGNNCTGTAGXGTCAGTGCTACGTACTNN\n>\nAGTTCTNNTAGTGCNNTNXNAGCCA\n";
    //                      ^^       ^^      ^                ^^           ^^      ^^ ^^^
    string color_data = "0\n1\n";

    vector<string> correct_out_seqs = {"AGTATTG", "CTGTAG", "GTCAGTGCTACGTACT", "AGTTCT", "TAGTGC", "T", "AGCCA"};
    vector<LL> correct_out_colors = {0,0,0,1,1,1,1};

    // Write to files

    string fasta1 = temp_file_manager.get_temp_file_name("fasta1");
    string colors1 = temp_file_manager.get_temp_file_name("colors1");
    
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

    cerr << correct_out_seqs << endl << out_seqs << endl << correct_out_colors << endl << out_colors << endl;
    assert(correct_out_seqs == out_seqs);
    assert(correct_out_colors == out_colors);

    cerr << "Delete non-ACGT test passed" << endl;

    return true;
}

bool misc_test_string_to_integer_safe(){
    try{
        string_to_integer_safe("1 2 3");
        assert(false); // Should have thrown an exception
    } catch(...){}
    try{
        string_to_integer_safe("  12 33 ");
        assert(false); // Should have thrown an exception
    } catch(...){}
    try{
        string_to_integer_safe("  12X33 ");
        assert(false); // Should have thrown an exception
    } catch(...){}
    try{
        string_to_integer_safe("  123XX ");
        assert(false); // Should have thrown an exception
    } catch(...){}
    try{
        string_to_integer_safe("   ");
        assert(false); // Should have thrown an exception
    } catch(...){}

    assert(string_to_integer_safe("  \n\t  \r 1234567890\n  \r\n") == 1234567890);

    cerr << "string_to_integer_safe test passed" << endl;
    return true;
}

bool misc_tests(){
    bool success = true;
    success &= misc_test_delete_non_ACGT();
    success &= misc_test_string_to_integer_safe();
    return success;
}