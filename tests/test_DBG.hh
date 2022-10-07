#pragma once

#include "setup_tests.hh"
#include <gtest/gtest.h>
#include "test_tools.hh"
#include "DBG.hh"
#include <unordered_set>

class DBG_Reference_Implementation{

public:

    vector<string> colex_kmers;
    unordered_map<string, vector<string>> outedges; // Outedges from a node are colex-sorted by destination
    unordered_map<string, vector<string>> inedges; // Inedges to a node are colex-sorted by source

    DBG_Reference_Implementation(){}
    DBG_Reference_Implementation(vector<string> reads, LL k){
        // Get sorted k-mers
        for(string S : reads) 

        std::sort(colex_kmers.begin(), colex_kmers.end(), colex_compare);
        colex_kmers.erase(unique(colex_kmers.begin(), colex_kmers.end()), colex_kmers.end()); // Erase duplicaes
        cout << colex_kmers << endl;

        std::unordered_set<std::string> kmers(colex_kmers.begin(), colex_kmers.end());

        // Get edges
        for(string S : reads){
            for(string x : get_all_distinct_kmers(S,k)){
                for(char c : string("ACGT")){
                    if(kmers.count(x.substr(1) + c))
                        outedges[x].push_back(x.substr(1) + c);
                    if(kmers.count(c + x.substr(0,k-1)))
                        inedges[x].push_back(c + x.substr(0,k-1));
                }
            }
        }
    }


};

class TEST_DBG : public ::testing::Test {
    protected:

    vector<string> reads;
    LL k;
    DBG_Reference_Implementation ref;
    plain_matrix_sbwt_t SBWT;
    SBWT_backward_traversal_support backward_support;
    DBG dbg;

    void SetUp() override {
        reads = {"CTGCGTAGTCGTACGATAAATTTCGATGTAGGCTCGTTCGGTCGC", "GACTTCTTTTCTTAGGCTAAAAAAAAA"};
        k = 3;
        ref = DBG_Reference_Implementation(reads,k);
        NodeBOSSInMemoryConstructor<plain_matrix_sbwt_t> builder;
        builder.build(reads, SBWT, k, true);
        backward_support = SBWT_backward_traversal_support(&SBWT);
        dbg = DBG(&SBWT);
    }

};

TEST_F(TEST_DBG, locate){
    // Check existing k-mers
    for(string kmer : ref.colex_kmers){
        DBG::Node v = dbg.locate(kmer);
        ASSERT_GE(v.id, 0); // Must be found
        ASSERT_EQ(dbg.get_node_label(v), kmer);
    }

    // check non-existent k-mer
    ASSERT_EQ(dbg.locate("CCG").id, -1);
}

TEST_F(TEST_DBG, iterate_all_nodes){
    LL kmer_idx = 0;
    for(DBG::Node v : dbg.all_nodes()){
        string fetched = backward_support.get_node_label(v.id);
        string correct = ref.colex_kmers[kmer_idx++];
        ASSERT_EQ(fetched, correct);
    }
    ASSERT_EQ(kmer_idx, ref.colex_kmers.size());
}

TEST_F(TEST_DBG, indegree){
    for(string kmer : ref.colex_kmers){
        DBG::Node v = dbg.locate(kmer);
        ASSERT_GE(v.id, 0); // Must be found
        if(dbg.indegree(v) != ref.inedges[kmer].size()){
            cout << "Break" << endl;
        }
        ASSERT_EQ(dbg.indegree(v), ref.inedges[kmer].size());
    }
}

TEST_F(TEST_DBG, inedges){
    for(DBG::Node v : dbg.all_nodes()){
        LL in_idx = 0;
        for(DBG::Edge e : dbg.inedges(v)){
            ASSERT_EQ(e.dest, v.id);
            string kmer_from = backward_support.get_node_label(e.source);
            string kmer_to = backward_support.get_node_label(e.dest);
            ASSERT_EQ(kmer_to.back(), e.label);
            ASSERT_EQ(kmer_from, ref.inedges[kmer_to][in_idx++]);
        }
    }
}

TEST_F(TEST_DBG, outdegree){
    for(string kmer : ref.colex_kmers){
        DBG::Node v = dbg.locate(kmer);
        ASSERT_GE(v.id, 0); // Must be found
        ASSERT_EQ(dbg.outdegree(v), ref.outedges[kmer].size());
    }
}


TEST_F(TEST_DBG, outedges){
    for(DBG::Node v : dbg.all_nodes()){
        LL out_idx = 0;
        for(DBG::Edge e : dbg.outedges(v)){
            ASSERT_EQ(e.source, v.id);
            string kmer_from = backward_support.get_node_label(e.source);
            string kmer_to = backward_support.get_node_label(e.dest);
            ASSERT_EQ(kmer_to.back(), e.label);
            ASSERT_EQ(kmer_to, ref.outedges[kmer_from][out_idx++]);
        }
    }
}