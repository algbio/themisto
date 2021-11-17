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

TEST(MISC_TEST, randomize_non_ACGT){
    vector<string> seqs = {"AATgCaTGCPPOjdjpqFbCL", "AACGTAAGCGALKJGF"};
    string fastafile = get_temp_file_manager().create_filename();
    write_as_fasta(seqs, fastafile);
    string fixedfile = fix_alphabet(fastafile, FASTA_MODE);
    Sequence_Reader_Buffered sr(fixedfile, FASTA_MODE);
    LL seqs_read = 0;

    auto is_dna = [](char c){
        c = toupper(c);
        return c == 'A' || c == 'C' || c == 'G' || c == 'T';
    };

    while(true){
        string S_fixed = sr.get_next_read();
        if(S_fixed.size() == 0) break;
        string S = seqs[seqs_read++];
        cout << S << endl << S_fixed << endl;
        ASSERT_EQ(S.size(), S_fixed.size());
        for(LL i = 0; i < S.size(); i++){
            if(is_dna(S[i])) {ASSERT_EQ(toupper(S[i]), toupper(S_fixed[i]));}
            if(!is_dna(S[i])) {ASSERT_TRUE(is_dna(S_fixed[i]));}
        }
    }
    ASSERT_EQ(seqs_read, seqs.size());
}

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
    Sequence_Reader_Buffered sr(fasta2, FASTA_MODE);
    while(true){
        string read = sr.get_next_read();
        if(read.size() > 0) out_seqs.push_back(read);
        else break;
    }

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
