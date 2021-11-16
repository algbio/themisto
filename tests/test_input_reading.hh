#pragma once

#include <gtest/gtest.h>
#include "../globals.hh"
#include "setup_tests.hh"

TEST(INPUT_PARSING, fasta_basic){
    vector<string> seqs = {"AAGTGCTGTANAYA","ACGTURYKMSWBDHVN-"};
    string fasta = ">\n" + seqs[0] + "\n>\n" + seqs[1] + "\n";
    string filename = string_to_temp_file(fasta);
    Sequence_Reader sr(filename, FASTA_MODE);
    ASSERT_FALSE(sr.done());
    ASSERT_EQ(sr.get_next_query_stream().get_all(), seqs[0]);
    ASSERT_FALSE(sr.done());
    ASSERT_EQ(sr.get_next_query_stream().get_all(), seqs[0]);
    ASSERT_TRUE(sr.done());

}

    // Check upper-casing

    // Check multiple lines

    // Check super long line (should not overflow any buffer)

//TEST(INPUT_PARSING, fastq){
    // Check upper-casing

    // Check multiple lines

    // Check super long line (should not overflow any buffer)

    // @-symbols in quality line

    // +-symbols in quality line


//}