#pragma once

#include <gtest/gtest.h>
#include "sbwt/stdlib_printing.hh"
#include "sbwt/SeqIO.hh"
#include "globals.hh"
#include "sbwt/globals.hh"
#include "setup_tests.hh"
#include "test_tools.hh"
#include "commands.hh"
#include "sbwt/variants.hh"
#include "sbwt/SBWT.hh"
#include "coloring/Coloring.hh"
#include "pseudoalign.hh"

using namespace sbwt;

struct CLI_test_files{
    string fastafile, indexprefix, tempdir, colorfile;
    CLI_test_files(vector<string> seqs, vector<int64_t> colors){
        fastafile = get_temp_file_manager().create_filename("",".fna");
        indexprefix = get_temp_file_manager().create_filename();
        tempdir = get_temp_file_manager().get_dir();
        colorfile = get_temp_file_manager().create_filename();

        write_as_fasta(seqs, fastafile);

        throwing_ofstream colors_out(colorfile);
        for(int64_t c : colors) colors_out << c << "\n";
    }
};


class CLI_TEST : public ::testing::Test {
    public:

    vector<string> seqs;
    vector<int64_t> colors;
    int64_t k;
    string fastafile, indexprefix, tempdir, colorfile;

    void SetUp() override {
        seqs = {"AACCGGTTA", "ACGTACGTG", "ATATGACATG"};
        colors = {3,1,2};
        k = 6;
        CLI_test_files f(seqs,colors);
        fastafile = f.fastafile;
        indexprefix = f.indexprefix;
        tempdir = f.tempdir;
        colorfile = f.colorfile;
    }

};

void load_sbwt_and_coloring(plain_matrix_sbwt_t& SBWT, Coloring<>& coloring, string indexprefix){
    SBWT.load(indexprefix + ".tdbg");
    coloring.load(indexprefix + ".tcolors", SBWT);
}

// If no colorfile is given, should assign colors automatically for each sequence
TEST_F(CLI_TEST, auto_colors){
    vector<string> args = {"build", "-k", to_string(k), "-i", fastafile, "-o", indexprefix, "--temp-dir", tempdir};
    cout << args << endl;
    sbwt::Argv argv(args);
    build_index_main(argv.size, argv.array);
    plain_matrix_sbwt_t SBWT; Coloring<> coloring;
    load_sbwt_and_coloring(SBWT, coloring, indexprefix);
    for(int64_t seq_id = 0; seq_id < seqs.size(); seq_id++){
        for(string kmer : get_all_kmers(seqs[seq_id], k)){
            int64_t node = SBWT.search(kmer);
            vector<int64_t> colorset = coloring.get_color_set_of_node_as_vector(node);
            ASSERT_EQ(colorset.size(), 1);
            ASSERT_EQ(colorset[0], seq_id);
        }
    }
}

TEST_F(CLI_TEST, reverse_complement_construction_with_auto_colors){
    vector<string> args = {"build", "-k", to_string(k), "-i", fastafile, "-o", indexprefix, "--temp-dir", tempdir, "--reverse-complements"};
    cout << args << endl;
    sbwt::Argv argv(args);
    build_index_main(argv.size, argv.array);
    plain_matrix_sbwt_t SBWT; Coloring<> coloring;
    load_sbwt_and_coloring(SBWT, coloring, indexprefix);
    for(int64_t seq_id = 0; seq_id < seqs.size(); seq_id++){
        for(string kmer : get_all_kmers(seqs[seq_id], k)){
            int64_t node = SBWT.search(kmer);
            ASSERT_GE(node, 0);
            vector<int64_t> colorset = coloring.get_color_set_of_node_as_vector(node);
            ASSERT_EQ(colorset.size(), 1);
            ASSERT_EQ(colorset[0], seq_id);
        }
        
        string rc = get_reverse_complement(seqs[seq_id]);
        for(string kmer : get_all_kmers(rc, k)){
            int64_t node = SBWT.search(kmer);
            ASSERT_GE(node, 0);
            vector<int64_t> colorset = coloring.get_color_set_of_node_as_vector(node);
            ASSERT_EQ(colorset.size(), 1);
            ASSERT_EQ(colorset[0], seq_id);
        }
    }
}

TEST_F(CLI_TEST, reverse_complement_construction_with_colorfile){
    vector<string> args = {"build", "-k", to_string(k), "-i", fastafile, "-o", indexprefix, "--temp-dir", tempdir, "--reverse-complements", "--color-file", colorfile};
    cout << args << endl;
    sbwt::Argv argv(args);
    build_index_main(argv.size, argv.array);
    plain_matrix_sbwt_t SBWT; Coloring<> coloring;
    load_sbwt_and_coloring(SBWT, coloring, indexprefix);
    for(int64_t seq_id = 0; seq_id < seqs.size(); seq_id++){
        for(string kmer : get_all_kmers(seqs[seq_id], k)){
            int64_t node = SBWT.search(kmer);
            ASSERT_GE(node, 0);
            vector<int64_t> colorset = coloring.get_color_set_of_node_as_vector(node);
            ASSERT_EQ(colorset.size(), 1);
            ASSERT_EQ(colorset[0], colors[seq_id]);
        }
        
        string rc = get_reverse_complement(seqs[seq_id]);
        for(string kmer : get_all_kmers(rc, k)){
            int64_t node = SBWT.search(kmer);
            ASSERT_GE(node, 0);
            vector<int64_t> colorset = coloring.get_color_set_of_node_as_vector(node);
            ASSERT_EQ(colorset.size(), 1);
            ASSERT_EQ(colorset[0], colors[seq_id]);
        }
    }
}

// If --no-colors is given, should not build colors. Also check that giving both --auto-colors and --color-file throws.
TEST_F(CLI_TEST, no_colors){
    vector<string> args = {"build", "--no-colors", "-k", to_string(k), "-i", fastafile, "-o", indexprefix, "--temp-dir", tempdir};
    sbwt::Argv argv(args);
    build_index_main(argv.size, argv.array);
    plain_matrix_sbwt_t SBWT; Coloring<> coloring;
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
    plain_matrix_sbwt_t SBWT; Coloring<> coloring;
    load_sbwt_and_coloring(SBWT, coloring, indexprefix);
    for(int64_t seq_id = 0; seq_id < seqs.size(); seq_id++){
        for(string kmer : get_all_kmers(seqs[seq_id], k)){
            int64_t node = SBWT.search(kmer);
            vector<int64_t> colorset = coloring.get_color_set_of_node_as_vector(node);
            ASSERT_EQ(colorset.size(), 1);
            ASSERT_EQ(colorset[0], colors[seq_id]);
        }
    }
}

TEST_F(CLI_TEST, gzip_input_in_building){
    // Create a gzipped file

    string gzip_outfile = get_temp_file_manager().create_filename("", ".fna.gz");
    { // Artifical scope to make the gz writer go out of scope to flush the stream. The flush method in the class does not actually flush
        sbwt::SeqIO::Writer<zstr::ofstream> gzip_out(gzip_outfile);
        for(string seq : seqs){
            gzip_out.write_sequence(seq.c_str(), seq.size());
        }
    }

    vector<string> args = {"build", "-k", to_string(k), "-i", gzip_outfile, "-o", indexprefix, "--temp-dir", tempdir};
    sbwt::Argv argv(args);
    build_index_main(argv.size, argv.array);
    plain_matrix_sbwt_t SBWT; Coloring<> coloring;
    load_sbwt_and_coloring(SBWT, coloring, indexprefix);
    for(int64_t seq_id = 0; seq_id < seqs.size(); seq_id++){
        for(string kmer : get_all_kmers(seqs[seq_id], k)){
            int64_t node = SBWT.search(kmer);
            vector<int64_t> colorset = coloring.get_color_set_of_node_as_vector(node);
            ASSERT_EQ(colorset.size(), 1);
            ASSERT_EQ(colorset[0], seq_id);
        }
    }
}

TEST(PREPROCESSING, upper_case){
    int64_t k = 4;

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

    plain_matrix_sbwt_t SBWT_1; Coloring<> coloring_1;
    plain_matrix_sbwt_t SBWT_2; Coloring<> coloring_2;

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

    vector<vector<int64_t> > res1 = parse_pseudoalignment_output_format_from_disk(q_resultfile1);
    vector<vector<int64_t> > res2 = parse_pseudoalignment_output_format_from_disk(q_resultfile2);

    ASSERT_EQ(res1, res2);
    vector<vector<int64_t> > correct = {{0},{0}};
    ASSERT_EQ(res1, correct);


}

TEST_F(CLI_TEST, test_color_matrix_dump){
    vector<string> args = {"build", "-k", to_string(k), "-i", fastafile, "-o", indexprefix, "--temp-dir", tempdir, "--color-file", colorfile};
    cout << args << endl;
    sbwt::Argv argv(args);
    build_index_main(argv.size, argv.array);
    plain_matrix_sbwt_t SBWT; Coloring<> coloring;
    load_sbwt_and_coloring(SBWT, coloring, indexprefix);

    // Do a dump in sparse format
    string sparse_dump_file = get_temp_file_manager().create_filename("", ".txt");
    vector<string> args2 = {"dump-color-matrix", "-i", indexprefix, "-o", sparse_dump_file, "--sparse"};
    sbwt::Argv argv2(args2);
    dump_color_matrix_main(argv2.size, argv2.array);

    // Do a dump in dense format
    string dense_dump_file = get_temp_file_manager().create_filename("", ".txt");
    vector<string> args3 = {"dump-color-matrix", "-i", indexprefix, "-o", dense_dump_file};
    sbwt::Argv argv3(args3);
    dump_color_matrix_main(argv3.size, argv3.array);

    // Extract true k-mers from the DBG (assuming the DBG works correctly)
    DBG dbg(&SBWT);
    vector<string> all_kmers;
    for(DBG::Node v : dbg.all_nodes()) 
        all_kmers.push_back(dbg.get_node_label(v));

    // Build true colors
    map<string, vector<int64_t> > true_colors; // kmer -> color set
    for(int64_t seq_idx = 0; seq_idx < seqs.size(); seq_idx++){
        string S = seqs[seq_idx];
        int64_t color = colors[seq_idx];
        for(string x : get_all_distinct_kmers(S, k)){
            true_colors[x].push_back(color);
        }
    }

    // Check the sparse dump
    string line;
    throwing_ifstream sparse_dump_in(sparse_dump_file);
    int64_t line_idx = 0;
    while(getline(sparse_dump_in.stream, line)){
        logger << "Line: " << line << endl;
        vector<string> tokens = split(line);
        string kmer = tokens[0];
        ASSERT_EQ(kmer, all_kmers[line_idx]); // Check that the k-mer is correct

        vector<int64_t> parsed_colors;
        for(int64_t i = 1; i < tokens.size(); i++){
            parsed_colors.push_back(stoll(tokens[i]));
        }

        ASSERT_EQ(true_colors[kmer], parsed_colors); // Check that the colors are correct
        line_idx++;
    }

    // Check the dense dump
    throwing_ifstream dense_dump_in(dense_dump_file);
    line_idx = 0;
    while(getline(dense_dump_in.stream, line)){
        logger << "Line: " << line << endl;
        vector<string> tokens = split(line);
        ASSERT_EQ(tokens.size(), 2);
        string kmer = tokens[0];
        ASSERT_EQ(kmer, all_kmers[line_idx]); // Check that the k-mer is correct

        vector<int64_t> parsed_colors;
        string ASCII_row = tokens[1];
        for(int64_t i = 0; i < ASCII_row.size(); i++){
            if(ASCII_row[i] == '1') parsed_colors.push_back(i);
        }

        ASSERT_EQ(true_colors[kmer], parsed_colors); // Check that the colors are correct
        line_idx++;
    }
}

TEST_F(CLI_TEST, multiple_input_files_user_colors){
    string gzip_outfile = get_temp_file_manager().create_filename("", ".fna.gz");

    int64_t m = 20; // Number of sequences
    int64_t k = 6;

    // Generate data
    vector<string> seqs;
    vector<int64_t> colors;
    for(int64_t i = 0; i < m; i++){
        seqs.push_back(get_random_dna_string(10,4));
        colors.push_back(i / 3); // Will be equal to file colors
    }

    // Split to files, 3 sequences per file
    vector<string> seqfiles;
    vector<string> colorfiles;
    for(int64_t i = 0; i < m; i += 3){
        vector<string> S(seqs.begin() + i, seqs.begin() + min(i + 3, m));
        vector<int64_t> C(colors.begin() + i, colors.begin() + min(i + 3, m));

        seqfiles.push_back(get_temp_file_manager().create_filename("", ".fna"));
        colorfiles.push_back(get_temp_file_manager().create_filename("", ".txt"));

        write_as_fasta(S, seqfiles.back());
        write_vector(C, colorfiles.back());
    }

    string seqfile_list = get_temp_file_manager().create_filename("", ".txt");
    string colorfile_list = get_temp_file_manager().create_filename("", ".txt");

    write_vector(seqfiles, seqfile_list);
    write_vector(colorfiles, colorfile_list);

    string index_prefix = get_temp_file_manager().create_filename();

    // Build from multiple files
    vector<string> args = {"build", "-k", to_string(k), "-i", seqfile_list, "-c", colorfile_list, "-o", index_prefix, "--temp-dir", tempdir};
    sbwt::Argv argv(args);
    build_index_main(argv.size, argv.array);

    // Build concatenated file
    string all_seq_file = get_temp_file_manager().create_filename("", ".fna");
    string all_colors_file = get_temp_file_manager().create_filename("", ".txt");
    write_as_fasta(seqs, all_seq_file);
    write_vector(colors, all_colors_file);

    string index_prefix2 = get_temp_file_manager().create_filename();
    vector<string> args2 = {"build", "-k", to_string(k), "-i", all_seq_file, "-c", all_colors_file, "-o", index_prefix2, "--temp-dir", tempdir};
    sbwt::Argv argv2(args2);
    build_index_main(argv2.size, argv2.array);

    // Build with file colors
    string index_prefix3 = get_temp_file_manager().create_filename();
    vector<string> args3 = {"build", "-k", to_string(k), "-i", seqfile_list, "--file-colors", "-o", index_prefix3, "--temp-dir", tempdir};
    sbwt::Argv argv3(args3);
    build_index_main(argv3.size, argv3.array);

    ASSERT_TRUE(files_are_equal(index_prefix + ".tdbg", index_prefix2 + ".tdbg"));
    ASSERT_TRUE(files_are_equal(index_prefix + ".tcolors", index_prefix2 + ".tcolors"));
    ASSERT_TRUE(files_are_equal(index_prefix2 + ".tcolors", index_prefix3 + ".tcolors"));

    // Check that file colors gives the same answer
}