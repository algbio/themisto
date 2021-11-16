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

// If no colorfile is given, should assign colors automatically for each sequence
TEST(MISC_TEST, cli_auto_colors){
    vector<string> seqs = {"AACCGGTT", "ACGTACGT", "ATATATAT"};
    LL k = 3;
    string fastafile = get_temp_file_manager().create_filename("",".fna");
    string indexdir = get_temp_file_manager().create_filename();
    string tempdir = get_temp_file_manager().get_dir();
    write_as_fasta(seqs, fastafile);
    vector<string> args = {"build", "-k", to_string(k), "-i", fastafile, "-o", indexdir, "--temp-dir", tempdir};
    Argv argv(args);
    build_index_main(argv.size, argv.array);
    Themisto themisto;
    themisto.load_from_directory(indexdir);
    for(LL seq_id = 0; seq_id < seqs.size(); seq_id++){
        for(string kmer : get_all_kmers(seqs[seq_id], k)){
            LL node = themisto.boss.find_kmer(kmer);
            vector<LL> colorset = themisto.coloring.get_colorvec(node, themisto.boss);
            ASSERT_EQ(colorset.size(), 1);
            ASSERT_EQ(colorset[0], seq_id);
        }
    }
}

// If --no-colors is given, should not build colors. Also check that giving both --auto-colors and --color-file throws.
TEST(MISC_TEST, no_colors){
    vector<string> seqs = {"AACCGGTT", "ACGTACGT", "ATATATAT"};
    LL k = 3;
    string fastafile = get_temp_file_manager().create_filename("",".fna");
    string indexdir = get_temp_file_manager().create_filename();
    string tempdir = get_temp_file_manager().get_dir();
    write_as_fasta(seqs, fastafile);
    vector<string> args = {"build", "--no-colors", "-k", to_string(k), "-i", fastafile, "-o", indexdir, "--temp-dir", tempdir};
    Argv argv(args);
    build_index_main(argv.size, argv.array);
    Themisto themisto;
    themisto.load_boss_from_directory(indexdir); // Should work
    try{
        themisto.load_colors_from_directory(indexdir); // Should throw
        FAIL(); // Did not throw
    } catch (const std::runtime_error &e){
        ASSERT_EQ(string(e.what()), "Error loading color data structure");
    }

    // Test --no-colors and --color-file at the same time

    string colorfile = get_temp_file_manager().create_filename();
    throwing_ofstream colors_out(colorfile);
    colors_out << "1\n2\n3\n";
    colors_out.close();
    vector<string> args_bad = {"build", "--no-colors", "--color-file", colorfile, "-k", to_string(k), "-i", fastafile, "-o", indexdir, "--temp-dir", tempdir};
    Argv argv_bad(args_bad);
    try{
        build_index_main(argv_bad.size, argv_bad.array);
        FAIL(); // Did not throw
    } catch(const std::runtime_error &e){
        ASSERT_EQ(string(e.what()), "Must not give both --no-colors and --colorfile");
    }
}