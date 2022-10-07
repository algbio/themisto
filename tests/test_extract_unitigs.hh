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
#include "new_coloring.hh"

using namespace sbwt;

// This is designed for k = 30
void construct_unitig_extraction_test_input(string fastafile, string colorfile){
    string random_data = "AATACCATGTCACAGCGTCAACGTTCAACACGCCCATACTGTGTTCCTGCGATGGGGGAGTAGTGACCCTCGATGCTGGACAATAACCCGATGACTAGACATAACGTCAATCTCGCTCTGTGTATTCGTATCGCCCATAGCTCTTCCATAAGCGTATTGAGTTGGGCTGTAAACACGGCCCCTTATAATTCTCTATTCTATGCACCGGTACCTATCTCAGGGCACAACTCCTGCCCGCTTTGGATTACTCGAGTACTGGCGCCACCTAATGCAGTACCCCCGAGTGGGATGAGTCAAATTTACGTGACCGGGAACACCATAGGTCCGCGTAAATACTGTGGGGCTATCCTTGGGCAACCTACTGTCACAGAGCTGGTACTCATATCTACATCACGCGTCGCAGAACAGCCAATCGGGCTGGATGTCAAAAGTAACAAGCGTGGTCCCTTAGGCGAAACCCGATCCTGATTTCAATAGGTTCCGTCCCGGGGGAGTAGCATCGGGCACGCAGCTTCCAAAACAGAATCGCCGTGGCCATATCGATGTCACGACGACGTTCACGGTTGGTTCCTCTTGAAGGCTCCCAGTCCTAGTGACGCGAGACGGTGATGTCGTGAGCGGAAGTGAACTCGGTCTTTGATTAATGTCAAAGCGCCGAGGCCCACGCATTCCCATCCACAAGTGTCTCTATGTGAGTGGTTTGTCCGCAAAGTAACCGGCGGACGCCCATCCCCGATCCTATAGCCGGAATAGTAGATGTTATAATTCTGAGGTATCGCCGTTGAGAGCTTATGCACCTCGGCCGTAGGGAGGTGGGAAGCGTTAGCGTCCCATAGACGGGGGTATATTTTCATCATGACTGTTGAAATCTGCGTTCGGAGGTTACACAGGGACGGAAGTGACTAATCGTGTCAAAGGACTTGTCTTCTCCTCTGCTCATGGGAAACCTACCACCACAGTCCCTGTGATGAAACGGATTTTCTCACTGTAGCTTTCAATT";

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

    write_as_fasta(seqs, fastafile);
    write_vector(colors, colorfile);

    logger << seqs << endl;

    ASSERT_EQ(seqs.size(), colors.size());
}

typedef uint32_t color_t;

class EXTRACT_UNITIGS_TEST : public testing::Test {

    private:

    vector<string> run_and_return_unitigs(string indexprefix, LL split_by_colors){
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
                for(LL i = 1; i < (LL)tokens.size(); i++){ // First token is unitig id -> skip
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
    Coloring coloring;
    sdsl::bit_vector is_dummy;
    vector<string> unitigs_with_colorsplit;
    vector<vector<color_t>> unitig_colors;
    vector<string> unitigs_without_colorsplit;
    DBG dbg;

    void SetUp() override {
        
        set_log_level(LogLevel::OFF);

        string seqfile = get_temp_file_manager().create_filename("",".fna");
        string colorfile = get_temp_file_manager().create_filename("",".txt");
        construct_unitig_extraction_test_input(seqfile, colorfile);
        LL k = 30;

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

        dbg = DBG(&SBWT);

        set_log_level(LogLevel::MAJOR);

    }

    virtual void TearDown() {}


};

// In all the tests below, when we write "indegree", it means the number of incoming
// edges from non-dummy nodes.

TEST_F(EXTRACT_UNITIGS_TEST, no_branches){
    LL k = SBWT.get_k();
    for(string unitig : unitigs_without_colorsplit){
        if(unitig.size() == k) continue; // Only one node: trivially ok
        for(LL i = 0; i < (LL)unitig.size()-k+1; i++){
            string kmer = unitig.substr(i,k);
            DBG::Node node = dbg.locate(kmer);
            ASSERT_GE(node.id, 0); // Exists in graph
            if(i == 0){ // Start node
                ASSERT_LT(dbg.outdegree(node), 2);
            }
            if(i > 0 && i < (LL)unitig.size()-k){ // Internal node
                ASSERT_LT(dbg.indegree(node), 2);
                ASSERT_LT(dbg.outdegree(node), 2);
            }
            if(i == (LL)unitig.size()-k){ // Last node
                ASSERT_LT(dbg.indegree(node), 2);
            }
            
        }
    }
}

TEST_F(EXTRACT_UNITIGS_TEST, split_by_colorsets){

    // Verify that the colorsets of all nodes in the unitig match the colorset of the unitig
    LL k = SBWT.get_k();
    ASSERT_EQ(unitigs_with_colorsplit.size(), unitig_colors.size());
    for(LL unitig_id = 0; unitig_id < (LL)unitigs_with_colorsplit.size(); unitig_id++){
        string unitig = unitigs_with_colorsplit[unitig_id];
        for(LL i = 0; i < (LL)unitig.size()-k+1; i++){
            LL node = SBWT.search(unitig.substr(i,k));
            vector<color_t> node_colors = coloring.get_color_set_as_vector(node);
            ASSERT_EQ(node_colors, unitig_colors[unitig_id]);
        }
    }
}

bool is_cyclic(DBG::Node first, DBG::Node last, DBG& dbg){
    if(dbg.indegree(first) != 1) return false;
    if(dbg.outdegree(last) != 1) return false;

    for(DBG::Edge e : dbg.inedges(first)) // There is exactly one in-edge
        return e.source == last.id;
    throw std::runtime_error("SHOULD NEVER COME HERE!!!!!");
}

TEST_F(EXTRACT_UNITIGS_TEST, maximality_no_color_split){
    LL k = SBWT.get_k();

    for(LL unitig_id = 0; unitig_id < (LL)unitigs_without_colorsplit.size(); unitig_id++){
        string unitig = unitigs_without_colorsplit[unitig_id];

        DBG::Node first = dbg.locate(unitig.substr(0, k));
        ASSERT_GE(first.id, 0); // Must exist in graph

        DBG::Node last = dbg.locate(unitig.substr((LL)unitig.size()-k, k));
        ASSERT_GE(last.id, 0); // Must exist in graph

        if(!is_cyclic(first, last, dbg)){
            // The indegree of the first node must be 0 or >= 2, or otherwise the predecessor must be forward-branching
            if(dbg.indegree(first) == 1){ // If indegree != 1, we are good
                DBG::Node pred = dbg.pred(first);
                ASSERT_TRUE(dbg.outdegree(pred) >= 2);
            }
            // The outdegree of the last node must be 0 or >= 2, or otherwise the successor must be backward-branching
            if(dbg.outdegree(last) == 1){ // If outdegree != 1, we are good
                DBG::Node succ = dbg.succ(last);
                ASSERT_TRUE(dbg.indegree(succ) >= 2);
            }
        }
    }
}

TEST_F(EXTRACT_UNITIGS_TEST, maximality_with_color_split){
    LL k = SBWT.get_k();

    for(LL unitig_id = 0; unitig_id < (LL)unitigs_with_colorsplit.size(); unitig_id++){
        string unitig = unitigs_with_colorsplit[unitig_id];

        DBG::Node first = dbg.locate(unitig.substr(0, k));
        ASSERT_GE(first.id, 0); // Must exist in graph

        DBG::Node last = dbg.locate(unitig.substr((LL)unitig.size()-k, k));
        ASSERT_GE(last.id, 0); // Must exist in graph

        if(!is_cyclic(first, last, dbg)){
            // The indegree of the first node must be 0 or >= 2, or otherwise the predecessor must be forward-branching
            if(dbg.indegree(first) == 1){ // If indegree != 1, we are good
                DBG::Node pred = dbg.pred(first);
                if(dbg.outdegree(pred) == 1){
                    vector<color_t> A = coloring.get_color_set_as_vector(first.id);
                    vector<color_t> B = coloring.get_color_set_as_vector(pred.id);
                    ASSERT_NE(A,B);
                }
            }
            // The outdegree of the last node must be 0 or >= 2, or otherwise the successor must be backward-branching
            if(dbg.outdegree(last) == 1){ // If outdegree != 1, we are good
                DBG::Node succ = dbg.succ(last);
                if(dbg.indegree(succ) == 1){
                    vector<color_t> A = coloring.get_color_set_as_vector(last.id);
                    vector<color_t> B = coloring.get_color_set_as_vector(succ.id);
                    ASSERT_NE(A,B);
                }
            }
        }
    }
}


// Check that every non-dummy node is in exactly one unitig
TEST_F(EXTRACT_UNITIGS_TEST, partition_without_colorsplit){
    LL k = SBWT.get_k();

    sdsl::bit_vector found = is_dummy; // dummies are marked as found from the beginning
    for(string unitig : unitigs_without_colorsplit){
        for(LL i = 0; i < (LL)unitig.size()-k+1; i++){
            DBG::Node node = dbg.locate(unitig.substr(i,k));
            ASSERT_EQ(found[node.id], 0);
            found[node.id] = 1;
        }
    }

    // Check that all were found
    for(LL i = 0; i < found.size(); i++){
        ASSERT_EQ(found[i], 1);
    }
}

// Check that every non-dummy node is in exactly one unitig
TEST_F(EXTRACT_UNITIGS_TEST, partition_with_colorsplit){
    LL k = SBWT.get_k();

    sdsl::bit_vector found = is_dummy; // dummies are marked as found from the beginning
    for(string unitig : unitigs_with_colorsplit){
        for(LL i = 0; i < (LL)unitig.size()-k+1; i++){
            DBG::Node node = dbg.locate(unitig.substr(i,k));
            ASSERT_EQ(found[node.id], 0);
            found[node.id] = 1;
        }
    }

    // Check that all were found
    for(LL i = 0; i < found.size(); i++){
        ASSERT_EQ(found[i], 1);
    }
}