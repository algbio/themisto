#pragma once

#include <gtest/gtest.h>
#include "sbwt/stdlib_printing.hh"
#include "globals.hh"
#include "sbwt/globals.hh"
#include "setup_tests.hh"
#include "test_tools.hh"
#include "commands.hh"
#include "sbwt/variants.hh"
#include "sbwt/SBWT.hh"
#include "new_coloring.hh"
#include "pseudoalign.hh"

using namespace sbwt;

typedef long long LL;

struct CLI_test_files{
    string fastafile, indexprefix, tempdir, colorfile;
    CLI_test_files(vector<string> seqs, vector<LL> colors){
        fastafile = get_temp_file_manager().create_filename("",".fna");
        indexprefix = get_temp_file_manager().create_filename();
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
    string fastafile, indexprefix, tempdir, colorfile;

    void SetUp() override {
        seqs = {"AACCGGTT", "ACGTACGT", "ATATATAT"};
        colors = {3,1,2};
        k = 3;
        CLI_test_files f(seqs,colors);
        fastafile = f.fastafile;
        indexprefix = f.indexprefix;
        tempdir = f.tempdir;
        colorfile = f.colorfile;
    }

};

void load_sbwt_and_coloring(plain_matrix_sbwt_t& SBWT, Coloring& coloring, string indexprefix){
    SBWT.load(indexprefix + ".tdbg");
    coloring.load(indexprefix + ".tcolors", SBWT);
}

// If no colorfile is given, should assign colors automatically for each sequence
TEST_F(CLI_TEST, auto_colors){
    vector<string> args = {"build", "-k", to_string(k), "-i", fastafile, "-o", indexprefix, "--temp-dir", tempdir};
    cout << args << endl;
    sbwt::Argv argv(args);
    build_index_main(argv.size, argv.array);
    plain_matrix_sbwt_t SBWT; Coloring coloring;
    load_sbwt_and_coloring(SBWT, coloring, indexprefix);
    for(LL seq_id = 0; seq_id < seqs.size(); seq_id++){
        for(string kmer : get_all_kmers(seqs[seq_id], k)){
            LL node = SBWT.search(kmer);
            vector<uint32_t> colorset = coloring.get_color_set_as_vector(node);
            ASSERT_EQ(colorset.size(), 1);
            ASSERT_EQ(colorset[0], seq_id);
        }
    }
}

// If --no-colors is given, should not build colors. Also check that giving both --auto-colors and --color-file throws.
TEST_F(CLI_TEST, no_colors){
    vector<string> args = {"build", "--no-colors", "-k", to_string(k), "-i", fastafile, "-o", indexprefix, "--temp-dir", tempdir};
    sbwt::Argv argv(args);
    build_index_main(argv.size, argv.array);
    plain_matrix_sbwt_t SBWT; Coloring coloring;
    SBWT.load(indexprefix + ".tdbg"); // Should work
    try{
        coloring.load(indexprefix, SBWT); // Should throw
        FAIL(); // Did not throw
    } catch (const std::exception &e){
        // Should come here
    }

    // Test --no-colors and --color-file at the same time
    vector<string> args_bad = {"build", "--no-colors", "--color-file", colorfile, "-k", to_string(k), "-i", fastafile, "-o", indexprefix, "--temp-dir", tempdir};
    sbwt::Argv argv_bad(args_bad);
    try{
        build_index_main(argv_bad.size, argv_bad.array);
        FAIL(); // Did not throw
    } catch(const std::exception &e){
        ASSERT_EQ(string(e.what()), "Must not give both --no-colors and --colorfile");
    }
}

TEST_F(CLI_TEST, build_colors_separately){

    // Build DBG
    vector<string> args = {"build", "--no-colors", "-k", to_string(k), "-i", fastafile, "-o", indexprefix, "--temp-dir", tempdir};
    sbwt::Argv argv(args);
    build_index_main(argv.size, argv.array);

    // Build Colors
    vector<string> args2 = {"build", "--load-dbg", "-c", colorfile, "-k", to_string(k), "-i", fastafile, "-o", indexprefix, "--temp-dir", tempdir};
    sbwt::Argv argv2(args2);
    build_index_main(argv2.size, argv2.array);

    // Check colors
    plain_matrix_sbwt_t SBWT; Coloring coloring;
    load_sbwt_and_coloring(SBWT, coloring, indexprefix);
    for(LL seq_id = 0; seq_id < seqs.size(); seq_id++){
        for(string kmer : get_all_kmers(seqs[seq_id], k)){
            LL node = SBWT.search(kmer);
            vector<uint32_t> colorset = coloring.get_color_set_as_vector(node);
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

    vector<string> args1 = {"build", "-k", to_string(k), "-i", f1.fastafile, "-c", f1.colorfile, "-o", f1.indexprefix, "--temp-dir", f1.tempdir};
    vector<string> args2 = {"build", "-k", to_string(k), "-i", f2.fastafile, "-c", f2.colorfile, "-o", f2.indexprefix, "--temp-dir", f2.tempdir};

    build_index_main(sbwt::Argv(args1).size, sbwt::Argv(args1).array);
    build_index_main(sbwt::Argv(args2).size, sbwt::Argv(args2).array);

    ASSERT_TRUE(files_are_equal(f1.indexprefix + ".tdbg", f2.indexprefix + ".tdbg"));
    ASSERT_TRUE(files_are_equal(f1.indexprefix + ".tcolors", f2.indexprefix + ".tcolors"));

    plain_matrix_sbwt_t SBWT_1; Coloring coloring_1;
    plain_matrix_sbwt_t SBWT_2; Coloring coloring_2;

    load_sbwt_and_coloring(SBWT_1, coloring_1, f1.indexprefix);
    load_sbwt_and_coloring(SBWT_2, coloring_2, f2.indexprefix);

    vector<string> queries = {seqs[0], seqs2[0]};

    string queryfile = get_temp_file_manager().create_filename("",".fna");
    write_as_fasta(queries, queryfile);

    string q_resultfile1 = get_temp_file_manager().create_filename();
    string q_resultfile2 = get_temp_file_manager().create_filename();

    vector<string> q_args1 = {"pseudoalign", "-q", queryfile, "-i", f1.indexprefix, "-o", q_resultfile1, "--temp-dir", f1.tempdir, "--rc"};
    vector<string> q_args2 = {"pseudoalign", "-q", queryfile, "-i", f2.indexprefix, "-o", q_resultfile2, "--temp-dir", f2.tempdir, "--rc"};

    pseudoalign_main(sbwt::Argv(q_args1).size, sbwt::Argv(q_args1).array);
    pseudoalign_main(sbwt::Argv(q_args2).size, sbwt::Argv(q_args2).array);

    ASSERT_TRUE(files_are_equal(q_resultfile1, q_resultfile2));

    vector<set<LL> > res1 = parse_pseudoalignment_output_format_from_disk(q_resultfile1);
    vector<set<LL> > res2 = parse_pseudoalignment_output_format_from_disk(q_resultfile2);

    ASSERT_EQ(res1, res2);
    vector<set<LL> > correct = {{0},{0}};
    ASSERT_EQ(res1, correct);
    

}