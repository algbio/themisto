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
#include "Themisto.hh"
#include "Argv.hh"
#include "extract_unitigs.hh"
#include <cassert>
#include <sstream>
#include "commands.hh"

class EXTRACT_UNITIGS_TEST_HAND_CRAFTED : public testing::Test {

    protected:

    string indexprefix;
    Themisto* themisto;
    vector<bool> is_dummy;
    vector<string> unitigs;

    void SetUp() override {
        string random_data = "AATACCATGTCACAGCGTCAACGTTCAACACGCCCATACTGTGTTCCTGCGATGGGGGAGTAGTGACCCTCGATGCTGGACAATAACCCGATGACTAGACATAACGTCAATCTCGCTCTGTGTATTCGTATCGCCCATAGCTCTTCCATAAGCGTATTGAGTTGGGCTGTAAACACGGCCCCTTATAATTCTCTATTCTATGCACCGGTACCTATCTCAGGGCACAACTCCTGCCCGCTTTGGATTACTCGAGTACTGGCGCCACCTAATGCAGTACCCCCGAGTGGGATGAGTCAAATTTACGTGACCGGGAACACCATAGGTCCGCGTAAATACTGTGGGGCTATCCTTGGGCAACCTACTGTCACAGAGCTGGTACTCATATCTACATCACGCGTCGCAGAACAGCCAATCGGGCTGGATGTCAAAAGTAACAAGCGTGGTCCCTTAGGCGAAACCCGATCCTGATTTCAATAGGTTCCGTCCCGGGGGAGTAGCATCGGGCACGCAGCTTCCAAAACAGAATCGCCGTGGCCATATCGATGTCACGACGACGTTCACGGTTGGTTCCTCTTGAAGGCTCCCAGTCCTAGTGACGCGAGACGGTGATGTCGTGAGCGGAAGTGAACTCGGTCTTTGATTAATGTCAAAGCGCCGAGGCCCACGCATTCCCATCCACAAGTGTCTCTATGTGAGTGGTTTGTCCGCAAAGTAACCGGCGGACGCCCATCCCCGATCCTATAGCCGGAATAGTAGATGTTATAATTCTGAGGTATCGCCGTTGAGAGCTTATGCACCTCGGCCGTAGGGAGGTGGGAAGCGTTAGCGTCCCATAGACGGGGGTATATTTTCATCATGACTGTTGAAATCTGCGTTCGGAGGTTACACAGGGACGGAAGTGACTAATCGTGTCAAAGGACTTGTCTTCTCCTCTGCTCATGGGAAACCTACCACCACAGTCCCTGTGATGAAACGGATTTTCTCACTGTAGCTTTCAATT";

        LL k = 30;

        vector<string> seqs;
        vector<LL> colors;

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

        ASSERT_EQ(seqs.size(), colors.size());

        themisto = new Themisto();
        string seqfile = get_temp_file_manager().create_filename("",".fna");
        string colorfile = get_temp_file_manager().create_filename("",".txt");
        write_as_fasta(seqs, seqfile);
        write_vector(colors, colorfile);
        
        set_log_level(LogLevel::MAJOR);

        // Build Themisto

        indexprefix = get_temp_file_manager().create_filename();

        stringstream argstring;
        argstring << "build -k"  << k << " --n-threads " << 4 << " --mem-megas " << 1024 << " -i " << seqfile << " -c " << colorfile << " -o " << indexprefix << " --temp-dir " << get_temp_file_manager().get_dir();
        Argv argv(split(argstring.str()));
        build_index_main(argv.size, argv.array);
        
        themisto->load(indexprefix);

        logger << "Getting dummy marks" << endl;
        is_dummy = themisto->boss.get_dummy_node_marks();
        
        // Call extract unitigs
        string unitigs_outfile = get_temp_file_manager().create_filename();
        stringstream argstring2;
        argstring2 << "extract-unitigs -i " << indexprefix << " --fasta-out " << unitigs_outfile;
        Argv argv2(split(argstring2.str()));
        extract_unitigs_main(argv2.size, argv2.array);

        throwing_ifstream unitigs_in(unitigs_outfile);
        // Parse unitigs
        string line;
        while(getline(unitigs_in.stream, line)){
            if(line[0] == '>') continue;
            unitigs.push_back(line);
        }

        set_log_level(LogLevel::OFF);

    }

    virtual void TearDown() {
        delete themisto;
    }


};

class EXTRACT_UNITIGS_TEST : public testing::Test {
   protected:

    static void SetUpTestCase() {
        
        // Avoid reallocating static objects if called in subclasses
        if (themisto == nullptr) {
            logger << "Building Themisto for extract_unitigs test suite" << endl;
            
            themisto = new Themisto();
            string seqfile = "example_input/coli3.fna";
            string colorfile = "example_input/colors.txt";
            LL k = 30;
            
            set_log_level(LogLevel::MAJOR);

            // Build Themisto

            indexprefix = get_temp_file_manager().create_filename();

            stringstream argstring;
            argstring << "build -k"  << k << " --n-threads " << 4 << " --mem-megas " << 1024 << " -i " << seqfile << " -c " << colorfile << " -o " << indexprefix << " --temp-dir " << get_temp_file_manager().get_dir();
            Argv argv(split(argstring.str()));
            build_index_main(argv.size, argv.array);
            
            themisto->load(indexprefix);

            logger << "Getting dummy marks" << endl;
            is_dummy = themisto->boss.get_dummy_node_marks();
            
            // Call extract unitigs
            string unitigs_outfile = get_temp_file_manager().create_filename();
            stringstream argstring2;
            argstring2 << "extract-unitigs -i " << indexprefix << " --fasta-out " << unitigs_outfile;
            Argv argv2(split(argstring2.str()));
            extract_unitigs_main(argv2.size, argv2.array);

            throwing_ifstream unitigs_in(unitigs_outfile);
            // Parse unitigs
            string line;
            while(getline(unitigs_in.stream, line)){
                if(line[0] == '>') continue;
                unitigs.push_back(line);
            }

            set_log_level(LogLevel::OFF);
        }
    }

    // Per-test-suite tear-down.
    // Called after the last test in this test suite.
    static void TearDownTestCase() {
        delete themisto;
        themisto = nullptr;
    }

    // Some expensive resources shared by all tests.
    static string indexprefix;
    static Themisto* themisto;
    static vector<string> unitigs;
    static vector<bool> is_dummy;
};

string EXTRACT_UNITIGS_TEST::indexprefix;
Themisto* EXTRACT_UNITIGS_TEST::themisto = nullptr;
vector<string> EXTRACT_UNITIGS_TEST::unitigs;
vector<bool> EXTRACT_UNITIGS_TEST::is_dummy;

// In all the tests below, when we write "indegree", it means the number of incoming
// edges from non-dummy nodes.


TEST_F(EXTRACT_UNITIGS_TEST_HAND_CRAFTED, no_branches){
    BOSS<sdsl::bit_vector>& boss = themisto->boss;
    LL k = boss.get_k();
    for(string unitig : unitigs){
        if(unitig.size() == k) continue; // Only one node: trivially ok
        for(LL i = 0; i < (LL)unitig.size()-k+1; i++){
            string kmer = unitig.substr(i,k);
            LL node = boss.find_kmer(kmer);
            ASSERT_GE(node, 0); // Exists in graph
            if(i == 0){ // Start node
                ASSERT_LT(boss.outdegree(node), 2);
            }
            if(i > 0 && i < (LL)unitig.size()-k){ // Internal node
                ASSERT_LT(boss.indegree(node), 2);
                ASSERT_LT(boss.outdegree(node), 2);
            }
            if(i == (LL)unitig.size()-k){ // Last node
                ASSERT_LT(boss.indegree(node), 2);
            }
            
        }
    }
}

TEST_F(EXTRACT_UNITIGS_TEST_HAND_CRAFTED, split_by_colorsets){
    // We're re-using the same Themisto index, so we use the same test fixture as
    // the other tests. However we don't case about the precomputed untigs.

    // Call extract with colorsets enabled
    string unitigs_outfile = get_temp_file_manager().create_filename();
    string colors_outfile = get_temp_file_manager().create_filename();
    stringstream argstring;
    argstring << "extract-unitigs -i " << indexprefix << " --fasta-out " << unitigs_outfile << " --colors-out " << colors_outfile;
    Argv argv(split(argstring.str()));
    extract_unitigs_main(argv.size, argv.array);

    vector<string> split_unitigs;
    vector<vector<LL> > split_unitig_colorsets;

    // Parse unitigs
    throwing_ifstream unitigs_in(unitigs_outfile);
    string line;
    while(getline(unitigs_in.stream, line)){
        if(line[0] == '>') continue;
        split_unitigs.push_back(line);
    }

    // Parse unitig colors
    throwing_ifstream colors_in(colors_outfile);
    while(getline(colors_in.stream, line)){
        vector<string> tokens = split(line, ' ');
        vector<LL> colors;
        for(LL i = 1; i < (LL)tokens.size(); i++){ // First token is unitig id -> skip
            colors.push_back(string_to_integer_safe(tokens[i]));
        }
        split_unitig_colorsets.push_back(colors);
    }

    // Verify that the colorsets of all nodes in the unitig match the colorset of the unitig
    LL k = themisto->boss.get_k();
    ASSERT_EQ(split_unitigs.size(), split_unitig_colorsets.size());
    for(LL unitig_id = 0; unitig_id < (LL)split_unitigs.size(); unitig_id++){
        ASSERT_GT(split_unitigs.size(), 0);
        string unitig = split_unitigs[unitig_id];
        for(LL i = 0; i < unitig.size()-k+1; i++){
            LL node = themisto->boss.find_kmer(unitig.substr(i,k));
            vector<LL> node_colors = themisto->coloring.get_colorvec(node, themisto->boss);
            ASSERT_EQ(node_colors, split_unitig_colorsets[unitig_id]);
        }
    }

    // todo for future: check that the split unitigs are maximal.
}

bool is_cyclic(LL first, LL last, BOSS<sdsl::bit_vector>& boss){
    if(boss.indegree(first) != 1) return false;
    if(boss.outdegree(last) != 1) return false;

    LL before_first = boss.edge_source(boss.inedge_range(first).first);
    if(before_first != last) return false;
    return true;
}

TEST_F(EXTRACT_UNITIGS_TEST_HAND_CRAFTED, maximality){
    BOSS<sdsl::bit_vector>& boss = themisto->boss;
    LL k = boss.get_k();
    for(LL unitig_id = 0; unitig_id < (LL)unitigs.size(); unitig_id++){
        string unitig = unitigs[unitig_id];

        LL first = boss.find_kmer(unitig.substr(0, k));
        ASSERT_GE(first, 0); // Must exist in graph

        LL last = boss.find_kmer(unitig.substr((LL)unitig.size()-k, k));
        ASSERT_GE(last, 0); // Must exist in graph

        if(!is_cyclic(first, last, boss)){
            // The indegree of the first node must be 0 or >= 2, or otherwise the predecessor must be forward-branching
            if(boss.indegree(first) == 1){ // If indegree != 1, we are good
                LL pred = boss.edge_source(boss.inedge_range(first).first); // predecessor
                ASSERT_TRUE(is_dummy[pred] || boss.outdegree(pred) >= 2);
            }
            // The outdegree of the last node must be 0 or >= 2, or otherwise the successor must be backward-branching
            if(boss.outdegree(last) == 1){ // If outdegree != 1, we are good
                LL succ = boss.walk(last, boss.outlabels_at(boss.outlabel_range(last).first));
                ASSERT_TRUE(boss.indegree(succ) >= 2);
            }
        }
    }
}

// Check that every non-dummy node is in exactly one unitig
TEST_F(EXTRACT_UNITIGS_TEST_HAND_CRAFTED, partition){
    BOSS<sdsl::bit_vector>& boss = themisto->boss;
    LL k = boss.get_k();
    vector<bool> found = is_dummy; // dummies are marked as found from the beginning
    for(string unitig : unitigs){
        for(LL i = 0; i < (LL)unitig.size()-k+1; i++){
            LL node = boss.find_kmer(unitig.substr(i,k));
            ASSERT_EQ(found[node], 0);
            found[node] = 1;
        }
    }

    // Check that all were found
    for(LL i = 0; i < found.size(); i++){
        ASSERT_EQ(found[i], 1);
    }

}
