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
            
            enable_logging();

            // Build Themisto

            string indexprefix = get_temp_file_manager().create_filename();

            stringstream argstring;
            argstring << "build -k"  << k << " --n-threads " << 4 << " --mem-megas " << 1024 << " -i " << seqfile << " -c " << colorfile << " -o " << indexprefix << " --temp-dir " << get_temp_file_manager().get_dir();
            Argv argv(split(argstring.str()));
            build_index_main(argv.size, argv.array);
            
            themisto->load(indexprefix);

            logger << "Getting dummy marks" << endl;
            is_dummy = themisto->boss.get_dummy_node_marks();
            
            // Extract unitigs and their colors

            stringstream unitigs_out;
            stringstream colors_out; // unused

            UnitigExtractor UE;
            UE.extract_unitigs(*themisto, unitigs_out, false, colors_out); // No color splits

            // Parse unitigs
            string line;
            while(getline(unitigs_out, line)){
                if(line[0] == '>') continue;
                unitigs.push_back(line);
            }

            disable_logging();
        }
    }

    // Per-test-suite tear-down.
    // Called after the last test in this test suite.
    static void TearDownTestCase() {
        delete themisto;
        themisto = nullptr;
    }

    // Some expensive resources shared by all tests.
    static Themisto* themisto;
    static vector<string> unitigs;
    static vector<bool> is_dummy;
};

Themisto* EXTRACT_UNITIGS_TEST::themisto = nullptr;
vector<string> EXTRACT_UNITIGS_TEST::unitigs;
vector<bool> EXTRACT_UNITIGS_TEST::is_dummy;

// In all the tests below, when we write "indegree", it means the number of incoming
// edges from non-dummy nodes.


TEST_F(EXTRACT_UNITIGS_TEST, no_branches){
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

TEST_F(EXTRACT_UNITIGS_TEST, split_by_colorsets){

    // First, compute the unitigs with color splits enabled

    BOSS<sdsl::bit_vector>& boss = themisto->boss;
    Coloring& coloring = themisto->coloring;
    UnitigExtractor UE;

    stringstream unitigs_out;
    stringstream colors_out;
    UE.extract_unitigs(*themisto, unitigs_out, true, colors_out);

    vector<string> split_unitigs;
    vector<vector<LL> > split_unitig_colorsets;

    // Parse unitigs
    string line;
    while(getline(unitigs_out, line)){
        if(line[0] == '>') continue;
        split_unitigs.push_back(line);
    }

    // Parse unitig colors
    while(getline(colors_out, line)){
        vector<string> tokens = split(line, ' ');
        vector<LL> colors;
        for(LL i = 1; i < (LL)tokens.size(); i++){ // First token is unitig id -> skip
            colors.push_back(string_to_integer_safe(tokens[i]));
        }
        split_unitig_colorsets.push_back(colors);
    }

    // Verify that the colorsets of all nodes in the unitig match the colorset of the unitig
    LL k = boss.get_k();
    ASSERT_EQ(split_unitigs.size(), split_unitig_colorsets.size());
    for(LL unitig_id = 0; unitig_id < (LL)split_unitigs.size(); unitig_id++){
        ASSERT_GT(split_unitigs.size(), 0);
        string unitig = split_unitigs[unitig_id];
        for(LL i = 0; i < unitig.size()-k+1; i++){
            LL node = boss.find_kmer(unitig.substr(i,k));
            vector<LL> node_colors = coloring.get_colorvec(node, boss);
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

TEST_F(EXTRACT_UNITIGS_TEST, maximality){
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
                LL succ = boss.walk(last, boss.outlabels_at(boss.outedge_range(last).first));
                ASSERT_TRUE(boss.indegree(succ) >= 2);
            }
        }
    }
}

// Check that every non-dummy node is in exactly one unitig
TEST_F(EXTRACT_UNITIGS_TEST, partition){
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
