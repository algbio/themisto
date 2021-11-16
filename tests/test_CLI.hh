#pragma once

#include <gtest/gtest.h>
#include "../stdlib_printing.hh"
#include "../globals.hh"
#include "setup_tests.hh"
#include "test_tools.hh"
#include "Themisto.hh"
#include "Argv.hh"
#include "commands.hh"

struct CLI_test_files{
    string fastafile, indexdir, tempdir, colorfile;
    CLI_test_files(vector<string> seqs, vector<LL> colors){
        fastafile = get_temp_file_manager().create_filename("",".fna");
        indexdir = get_temp_file_manager().create_filename();
        tempdir = get_temp_file_manager().get_dir();
        colorfile = get_temp_file_manager().create_filename();

        write_as_fasta(seqs, fastafile);

        throwing_ofstream colors_out(colorfile);
        for(LL c : colors) colors_out << c << "\n";
    }
};


class CLI_TEST : public ::testing::Test {
    public:

    vector<string> seqs;
    vector<LL> colors;
    LL k;
    string fastafile, indexdir, tempdir, colorfile;

    void SetUp() override {
        seqs = {"AACCGGTT", "ACGTACGT", "ATATATAT"};
        colors = {3,1,2};
        k = 3;
        CLI_test_files f(seqs,colors);
        fastafile = f.fastafile;
        indexdir = f.indexdir;
        tempdir = f.tempdir;
        colorfile = f.colorfile;
    }

};

// If no colorfile is given, should assign colors automatically for each sequence
TEST_F(CLI_TEST, auto_colors){
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
TEST_F(CLI_TEST, no_colors){
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
    vector<string> args_bad = {"build", "--no-colors", "--color-file", colorfile, "-k", to_string(k), "-i", fastafile, "-o", indexdir, "--temp-dir", tempdir};
    Argv argv_bad(args_bad);
    try{
        build_index_main(argv_bad.size, argv_bad.array);
        FAIL(); // Did not throw
    } catch(const std::runtime_error &e){
        ASSERT_EQ(string(e.what()), "Must not give both --no-colors and --colorfile");
    }
}

TEST_F(CLI_TEST, build_colors_separately){

    // Build DBG
    vector<string> args = {"build", "--no-colors", "-k", to_string(k), "-i", fastafile, "-o", indexdir, "--temp-dir", tempdir};
    Argv argv(args);
    build_index_main(argv.size, argv.array);

    // Build Colors
    vector<string> args2 = {"build", "--load-dbg", "-c", colorfile, "-k", to_string(k), "-i", fastafile, "-o", indexdir, "--temp-dir", tempdir};
    Argv argv2(args2);
    build_index_main(argv2.size, argv2.array);

    // Check colors
    Themisto themisto;
    themisto.load_from_directory(indexdir);
    for(LL seq_id = 0; seq_id < seqs.size(); seq_id++){
        for(string kmer : get_all_kmers(seqs[seq_id], k)){
            LL node = themisto.boss.find_kmer(kmer);
            vector<LL> colorset = themisto.coloring.get_colorvec(node, themisto.boss);
            ASSERT_EQ(colorset.size(), 1);
            ASSERT_EQ(colorset[0], colors[seq_id]);
        }
    }
}

TEST(PREPROCESSING, upper_case){
    LL k = 4;

    vector<string> seqs = {"AGGTCGATTCGATCGATGC"};
    vector<string> seqs2 = {"AGGtCGATTcGATCgaTGC"};

    CLI_test_files f1(seqs, {0});
    CLI_test_files f2(seqs2, {0});

    vector<string> args1 = {"build", "-k", to_string(k), "-i", f1.fastafile, "-c", f1.colorfile, "-o", f1.indexdir, "--temp-dir", f1.tempdir};
    vector<string> args2 = {"build", "-k", to_string(k), "-i", f2.fastafile, "-c", f2.colorfile, "-o", f2.indexdir, "--temp-dir", f2.tempdir};

    Argv argv1(args1);
    Argv argv2(args2);

    build_index_main(argv1.size, argv1.array);
    build_index_main(argv2.size, argv2.array);

    Themisto themisto1;
    Themisto themisto2;

    themisto1.load_from_directory(f1.indexdir);
    themisto2.load_from_directory(f2.indexdir);

    ASSERT_TRUE(themisto1.boss == themisto2.boss);

    vector<string> queries = {seqs[0], seqs2[0]};

    string queryfile = get_temp_file_manager().create_filename("",".fna");
    write_as_fasta(queries, queryfile);

    string q_resultfile1 = get_temp_file_manager().create_filename();
    string q_resultfile2 = get_temp_file_manager().create_filename();

    vector<string> q_args1 = {"pseudoalign", "-q", queryfile, "-i", f1.indexdir, "-o", q_resultfile1, "--temp-dir", f1.tempdir, "--rc"};
    vector<string> q_args2 = {"pseudoalign", "-q", queryfile, "-i", f2.indexdir, "-o", q_resultfile2, "--temp-dir", f2.tempdir, "--rc"};

    Argv q_argv1(q_args1);
    Argv q_argv2(q_args2);

    pseudoalign_main(q_argv1.size, q_argv1.array);
    pseudoalign_main(q_argv2.size, q_argv2.array);

    ASSERT_TRUE(files_are_equal(q_resultfile1, q_resultfile2));

    vector<set<LL> > res1 = Themisto::parse_output_format_from_disk(q_resultfile1);
    vector<set<LL> > res2 = Themisto::parse_output_format_from_disk(q_resultfile2);

    ASSERT_EQ(res1, res2);
    vector<set<LL> > correct = {{0},{0}};
    ASSERT_EQ(res1, correct);
    

}