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
#include "globals.hh"
#include "setup_tests.hh"
#include "sbwt/globals.hh"
#include "extract_unitigs.hh"
#include <cassert>
#include <sstream>
#include "commands.hh"
#include "DBG.hh"
#include "sbwt/SBWT.hh"
#include "coloring/Coloring.hh"
#include "coloring/Coloring_Builder.hh"

using namespace sbwt;

// This is designed for k = 30
void construct_unitig_extraction_test_input(string fastafile, string colorfile){
    string random_data = "AATACCATGTCACAGCGTCAACGTTCAACACGCCCATACTGTGTTCCTGCGATGGGGGAGTAGTGACCCTCGATGCTGGACAATAACCCGATGACTAGACATAACGTCAATCTCGCTCTGTGTATTCGTATCGCCCATAGCTCTTCCATAAGCGTATTGAGTTGGGCTGTAAACACGGCCCCTTATAATTCTCTATTCTATGCACCGGTACCTATCTCAGGGCACAACTCCTGCCCGCTTTGGATTACTCGAGTACTGGCGCCACCTAATGCAGTACCCCCGAGTGGGATGAGTCAAATTTACGTGACCGGGAACACCATAGGTCCGCGTAAATACTGTGGGGCTATCCTTGGGCAACCTACTGTCACAGAGCTGGTACTCATATCTACATCACGCGTCGCAGAACAGCCAATCGGGCTGGATGTCAAAAGTAACAAGCGTGGTCCCTTAGGCGAAACCCGATCCTGATTTCAATAGGTTCCGTCCCGGGGGAGTAGCATCGGGCACGCAGCTTCCAAAACAGAATCGCCGTGGCCATATCGATGTCACGACGACGTTCACGGTTGGTTCCTCTTGAAGGCTCCCAGTCCTAGTGACGCGAGACGGTGATGTCGTGAGCGGAAGTGAACTCGGTCTTTGATTAATGTCAAAGCGCCGAGGCCCACGCATTCCCATCCACAAGTGTCTCTATGTGAGTGGTTTGTCCGCAAAGTAACCGGCGGACGCCCATCCCCGATCCTATAGCCGGAATAGTAGATGTTATAATTCTGAGGTATCGCCGTTGAGAGCTTATGCACCTCGGCCGTAGGGAGGTGGGAAGCGTTAGCGTCCCATAGACGGGGGTATATTTTCATCATGACTGTTGAAATCTGCGTTCGGAGGTTACACAGGGACGGAAGTGACTAATCGTGTCAAAGGACTTGTCTTCTCCTCTGCTCATGGGAAACCTACCACCACAGTCCCTGTGATGAAACGGATTTTCTCACTGTAGCTTTCAATT";

    vector<string> seqs;
    vector<int64_t> colors;

    // Single nucleotide change

    string S = random_data.substr(0,150);
    seqs.push_back(S);
    colors.push_back(0);

    S[75] = (S[75] == 'A') ? 'C' : 'A';
    seqs.push_back(S);
    colors.push_back(0);

    // Self-loop
    seqs.push_back(random_data.substr(150,50)
                    + "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                    + random_data.substr(200,50));
    colors.push_back(0);

    // Long periodic string
    seqs.push_back(random_data.substr(250,20) + random_data.substr(250,20) + random_data.substr(250,20));
    colors.push_back(0);

    // Color ends in the middle of a unitig
    seqs.push_back(random_data.substr(270, 100));
    colors.push_back(1);
    seqs.push_back(random_data.substr(270, 80));
    colors.push_back(2);

    // Color starts in the middle of a unitig
    seqs.push_back(random_data.substr(370, 100));
    colors.push_back(3);
    seqs.push_back(random_data.substr(400, 100));
    colors.push_back(4);

    // Color fully contained in a unitig
    seqs.push_back(random_data.substr(500, 100));
    colors.push_back(5);
    seqs.push_back(random_data.substr(525, 50));
    colors.push_back(6);

    // Isolated unitig
    seqs.push_back(random_data.substr(600, 50));
    colors.push_back(7);

    // X-shape
    seqs.push_back("A" + random_data.substr(650, 29) + "AG");
    seqs.push_back("C" + random_data.substr(650, 29) + "AT");
    colors.push_back(0);
    colors.push_back(0);

    // Isolated k-mer
    seqs.push_back(random_data.substr(700,30));
    colors.push_back(0);

    // Shotert than k
    seqs.push_back(random_data.substr(730,20));
    colors.push_back(0);

    write_as_fasta(seqs, fastafile);
    write_vector(colors, colorfile);

    logger << seqs << endl;

    ASSERT_EQ(seqs.size(), colors.size());
}

typedef int64_t color_t;

class EXTRACT_UNITIGS_TEST : public testing::Test {

    private:

    vector<string> run_and_return_unitigs(string indexprefix, int64_t split_by_colors){
        // Call extract unitigs
        string unitigs_outfile = get_temp_file_manager().create_filename();
        string colors_outfile = get_temp_file_manager().create_filename();
        string gfa_outfile = get_temp_file_manager().create_filename();
        stringstream argstring2;
        argstring2 << "extract-unitigs -i " << indexprefix << " --fasta-out " << unitigs_outfile << " --gfa-out " << gfa_outfile;
        if(split_by_colors) argstring2 << " --colors-out " << colors_outfile;
        Argv argv2(split(argstring2.str()));
        extract_unitigs_main(argv2.size, argv2.array);

        // Parse unitigs
        throwing_ifstream unitigs_in(unitigs_outfile);

        string line;
        vector<string> parsed_unitigs;
        while(getline(unitigs_in.stream, line)){
            if(line[0] == '>') continue;
            parsed_unitigs.push_back(line);
        }

        if(split_by_colors){
            // Parse unitig colors
            throwing_ifstream colors_in(colors_outfile);
            while(getline(colors_in.stream, line)){
                vector<string> tokens = split(line, ' ');
                vector<color_t> colors;
                for(int64_t i = 1; i < (int64_t)tokens.size(); i++){ // First token is unitig id -> skip
                    colors.push_back(string_to_integer_safe(tokens[i]));
                }
                unitig_colors.push_back(colors);
            }
        }

        return parsed_unitigs;

    }

    protected:

    string indexprefix;
    plain_matrix_sbwt_t SBWT;
    Coloring<> coloring;
    sdsl::bit_vector is_dummy;
    vector<string> unitigs_with_colorsplit;
    vector<vector<color_t>> unitig_colors;
    vector<string> unitigs_without_colorsplit;
    DBG* dbg = nullptr;

    void SetUp() override {

        set_log_level(LogLevel::OFF);

        string seqfile = get_temp_file_manager().create_filename("",".fna");
        string colorfile = get_temp_file_manager().create_filename("",".txt");
        construct_unitig_extraction_test_input(seqfile, colorfile);
        int64_t k = 30;

        // Build Themisto

        indexprefix = get_temp_file_manager().create_filename();

        stringstream argstring;
        argstring << "build -k"  << k << " --n-threads " << 4 << " --mem-megas " << 2048 << " -i " << seqfile << " -c " << colorfile << " -o " << indexprefix << " --temp-dir " << get_temp_file_manager().get_dir();
        Argv argv(split(argstring.str()));
        build_index_main(argv.size, argv.array);

        // Compute unitigs, with and without the color split

        unitigs_without_colorsplit = run_and_return_unitigs(indexprefix, false);
        unitigs_with_colorsplit = run_and_return_unitigs(indexprefix, true);

        SBWT.load(indexprefix + ".tdbg");
        coloring.load(indexprefix + ".tcolors", SBWT);
        logger << "Getting dummy marks" << endl;
        is_dummy = SBWT.compute_dummy_node_marks();

        dbg = new DBG(&SBWT);

        set_log_level(LogLevel::MAJOR);

    }

    virtual void TearDown() {
        delete dbg;
    }


};

// In all the tests below, when we write "indegree", it means the number of incoming
// edges from non-dummy nodes.

TEST_F(EXTRACT_UNITIGS_TEST, no_branches){
    int64_t k = SBWT.get_k();
    for(string unitig : unitigs_without_colorsplit){
        if(unitig.size() == k) continue; // Only one node: trivially ok
        for(int64_t i = 0; i < (int64_t)unitig.size()-k+1; i++){
            string kmer = unitig.substr(i,k);
            DBG::Node node = dbg->locate(kmer);
            ASSERT_GE(node.id, 0); // Exists in graph
            if(i == 0){ // Start node
                ASSERT_LT(dbg->outdegree(node), 2);
            }
            if(i > 0 && i < (int64_t)unitig.size()-k){ // Internal node
                ASSERT_LT(dbg->indegree(node), 2);
                ASSERT_LT(dbg->outdegree(node), 2);
            }
            if(i == (int64_t)unitig.size()-k){ // Last node
                ASSERT_LT(dbg->indegree(node), 2);
            }

        }
    }
}

TEST_F(EXTRACT_UNITIGS_TEST, split_by_colorsets){

    // Verify that the colorsets of all nodes in the unitig match the colorset of the unitig
    int64_t k = SBWT.get_k();
    ASSERT_EQ(unitigs_with_colorsplit.size(), unitig_colors.size());
    for(int64_t unitig_id = 0; unitig_id < (int64_t)unitigs_with_colorsplit.size(); unitig_id++){
        string unitig = unitigs_with_colorsplit[unitig_id];
        for(int64_t i = 0; i < (int64_t)unitig.size()-k+1; i++){
            int64_t node = SBWT.search(unitig.substr(i,k));
            vector<color_t> node_colors = coloring.get_color_set_of_node_as_vector(node);
            ASSERT_EQ(node_colors, unitig_colors[unitig_id]);
        }
    }
}

bool is_cyclic(DBG::Node first, DBG::Node last, DBG* dbg){
    if(dbg->indegree(first) != 1) return false;
    if(dbg->outdegree(last) != 1) return false;

    for(DBG::Edge e : dbg->inedges(first)) // There is exactly one in-edge
        return e.source == last.id;
    throw std::runtime_error("SHOULD NEVER COME HERE!!!!!");
}

TEST_F(EXTRACT_UNITIGS_TEST, maximality_no_color_split){
    int64_t k = SBWT.get_k();

    for(int64_t unitig_id = 0; unitig_id < (int64_t)unitigs_without_colorsplit.size(); unitig_id++){
        string unitig = unitigs_without_colorsplit[unitig_id];

        DBG::Node first = dbg->locate(unitig.substr(0, k));
        ASSERT_GE(first.id, 0); // Must exist in graph

        DBG::Node last = dbg->locate(unitig.substr((int64_t)unitig.size()-k, k));
        ASSERT_GE(last.id, 0); // Must exist in graph

        if(!is_cyclic(first, last, dbg)){
            // The indegree of the first node must be 0 or >= 2, or otherwise the predecessor must be forward-branching
            if(dbg->indegree(first) == 1){ // If indegree != 1, we are good
                DBG::Node pred = dbg->pred(first);
                ASSERT_TRUE(dbg->outdegree(pred) >= 2);
            }
            // The outdegree of the last node must be 0 or >= 2, or otherwise the successor must be backward-branching
            if(dbg->outdegree(last) == 1){ // If outdegree != 1, we are good
                DBG::Node succ = dbg->succ(last);
                ASSERT_TRUE(dbg->indegree(succ) >= 2);
            }
        }
    }
}

TEST_F(EXTRACT_UNITIGS_TEST, maximality_with_color_split){
    int64_t k = SBWT.get_k();

    for(int64_t unitig_id = 0; unitig_id < (int64_t)unitigs_with_colorsplit.size(); unitig_id++){
        string unitig = unitigs_with_colorsplit[unitig_id];

        DBG::Node first = dbg->locate(unitig.substr(0, k));
        ASSERT_GE(first.id, 0); // Must exist in graph

        DBG::Node last = dbg->locate(unitig.substr((int64_t)unitig.size()-k, k));
        ASSERT_GE(last.id, 0); // Must exist in graph

        if(!is_cyclic(first, last, dbg)){
            // The indegree of the first node must be 0 or >= 2, or otherwise the predecessor must be forward-branching
            if(dbg->indegree(first) == 1){ // If indegree != 1, we are good
                DBG::Node pred = dbg->pred(first);
                if(dbg->outdegree(pred) == 1){
                    vector<color_t> A = coloring.get_color_set_of_node_as_vector(first.id);
                    vector<color_t> B = coloring.get_color_set_of_node_as_vector(pred.id);
                    ASSERT_NE(A,B);
                }
            }
            // The outdegree of the last node must be 0 or >= 2, or otherwise the successor must be backward-branching
            if(dbg->outdegree(last) == 1){ // If outdegree != 1, we are good
                DBG::Node succ = dbg->succ(last);
                if(dbg->indegree(succ) == 1){
                    vector<color_t> A = coloring.get_color_set_of_node_as_vector(last.id);
                    vector<color_t> B = coloring.get_color_set_of_node_as_vector(succ.id);
                    ASSERT_NE(A,B);
                }
            }
        }
    }
}


// Check that every non-dummy node is in exactly one unitig
TEST_F(EXTRACT_UNITIGS_TEST, partition_without_colorsplit){
    int64_t k = SBWT.get_k();

    sdsl::bit_vector found = is_dummy; // dummies are marked as found from the beginning
    for(string unitig : unitigs_without_colorsplit){
        for(int64_t i = 0; i < (int64_t)unitig.size()-k+1; i++){
            DBG::Node node = dbg->locate(unitig.substr(i,k));
            ASSERT_EQ(found[node.id], 0);
            found[node.id] = 1;
        }
    }

    // Check that all were found
    for(int64_t i = 0; i < found.size(); i++){
        ASSERT_EQ(found[i], 1);
    }
}

// Check that every non-dummy node is in exactly one unitig
TEST_F(EXTRACT_UNITIGS_TEST, partition_with_colorsplit){
    int64_t k = SBWT.get_k();

    sdsl::bit_vector found = is_dummy; // dummies are marked as found from the beginning
    for(string unitig : unitigs_with_colorsplit){
        for(int64_t i = 0; i < (int64_t)unitig.size()-k+1; i++){
            DBG::Node node = dbg->locate(unitig.substr(i,k));
            ASSERT_EQ(found[node.id], 0);
            found[node.id] = 1;
        }
    }

    // Check that all were found
    for(int64_t i = 0; i < found.size(); i++){
        ASSERT_EQ(found[i], 1);
    }
}

// Returns pair (unitigs, color sets)
pair<vector<string>, vector<vector<int64_t>>> get_colored_unitigs_with_themisto(string input_file_listfile, int64_t k){
    // Build Themisto with file colors
    string indexprefix = get_temp_file_manager().create_filename();
    vector<string> args = {"build", "-k", to_string(k), "-i", input_file_listfile, "-o", indexprefix, "--temp-dir", get_temp_file_manager().get_dir(), "--file-colors"};
    cout << args << endl;
    sbwt::Argv argv(args);
    build_index_main(argv.size, argv.array);
    plain_matrix_sbwt_t SBWT;
    SBWT.load(indexprefix + ".tdbg");
    Coloring<> coloring;
    coloring.load(indexprefix + ".tcolors", SBWT);

    UnitigExtractor<Coloring<SDSL_Variant_Color_Set>> UE;
    DBG dbg(&SBWT);

    string unitigs_outfile = get_temp_file_manager().create_filename("unitigs-",".fna");
    string unitig_colors_outfile = get_temp_file_manager().create_filename("unitigs-colors-",".txt");

    sbwt::throwing_ofstream unitigs_out(unitigs_outfile);
    sbwt::throwing_ofstream unitig_colors_out(unitig_colors_outfile);

    sbwt::SeqIO::NullStream gfa_null_stream;
    UE.extract_unitigs(dbg, coloring, unitigs_out.stream, true, unitig_colors_out.stream, gfa_null_stream, 0);

    unitigs_out.close();
    unitig_colors_out.close();

    // Parse unitigs and colors from disk
    vector<string> unitigs;
    vector<vector<int64_t> > color_sets;
    SeqIO::Reader<> unitigs_in(unitigs_outfile);
    while(true){ // Read unitigs
        string S = unitigs_in.get_next_read();
        if(S == "") break;
        else unitigs.push_back(S);
    }
    sbwt::throwing_ifstream color_sets_in(unitig_colors_outfile);
    string line;
    while(getline(color_sets_in.stream, line)){ // Read colors
        vector<string> tokens = split(line);
        vector<int64_t> colors;
        for(int64_t i = 1; i < tokens.size(); i++){ // i == 0 is the unitig id, which we ignore here
            colors.push_back(fast_string_to_int(tokens[i].c_str(), tokens[i].size()));
        }
        color_sets.push_back(colors);
    }

    return {unitigs, color_sets};
}

TEST(TEST_GGCAT, check_vs_themisto){
    string fastafile = get_temp_file_manager().create_filename("",".fna");
    string unused_colorfile = get_temp_file_manager().create_filename("",".txt");
    construct_unitig_extraction_test_input(fastafile, unused_colorfile); int64_t k = 30; // Putting k on the same line here because the test inputs are designed for k = 30

    // Split the unitigs into one unitig per file for ggcat
    vector<string> ggcat_input_files = split_seqs_to_separate_files(fastafile);
    string input_file_listfile = get_temp_file_manager().create_filename("", ".txt");
    write_lines(ggcat_input_files, input_file_listfile);

    // Run Themisto
    vector<string> themisto_unitigs;
    vector<vector<int64_t> > themisto_color_sets;
    std::tie(themisto_unitigs, themisto_color_sets) = get_colored_unitigs_with_themisto(input_file_listfile, k);

    // Run ggcat
    Colored_Unitig_Stream_GGCAT US_GGCAT(ggcat_input_files, 2, 3, k, false); // No reverse complements

    vector<string> ggcat_unitigs;
    vector<vector<int64_t>> ggcat_color_sets;
    while(!US_GGCAT.done()){
        ggcat_unitigs.push_back(US_GGCAT.next_unitig());
        ggcat_color_sets.push_back(US_GGCAT.next_colors());
    }

    ASSERT_EQ(ggcat_unitigs, themisto_unitigs);
    ASSERT_EQ(ggcat_color_sets, themisto_color_sets);


}