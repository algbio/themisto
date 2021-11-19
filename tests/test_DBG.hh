#pragma once

#include "setup_tests.hh"
#include <gtest/gtest.h>
#include "DBG.hh"
#include "BOSS_builder.hh"

class DBG_Reference_Implementation{

public:

    vector<string> colex_kmers;
    unordered_map<string, vector<string>> outedges; // Outedges from a node are colex-sorted by destination
    unordered_map<string, vector<string>> inedges; // Inedges to a node are colex-sorted by source

    DBG_Reference_Implementation(){}
    DBG_Reference_Implementation(vector<string> reads, LL k){
        // Get sorted k-mers
        for(string S : reads) 
            if(S.size() >= k+1) // Edge centric
                for(string x : get_all_distinct_kmers(S, k)) 
                    colex_kmers.push_back(x);
        std::sort(colex_kmers.begin(), colex_kmers.end(), colex_compare);
        colex_kmers.erase(unique(colex_kmers.begin(), colex_kmers.end()), colex_kmers.end()); // Erase duplicaes
        cout << colex_kmers << endl;

        // Get edges
        for(string S : reads){
            for(string x : get_all_distinct_kmers(S, k+1)){
                string source = x.substr(0,k);
                string dest = x.substr(1,k);
                outedges[source].push_back(dest);
                inedges[dest].push_back(source);
            }
        }

        // Sort out-edges and delete duplicates
        for(auto& keyval : outedges){
            vector<string>& out = keyval.second;
            std::sort(out.begin(), out.end(), colex_compare);
            out.erase(unique(out.begin(), out.end()), out.end()); // Erase duplicaes
        }

        // Sort in-edges and delete duplicates
        for(auto& keyval : inedges){
            vector<string>& in = keyval.second;
            std::sort(in.begin(), in.end(), colex_compare);
            in.erase(unique(in.begin(), in.end()), in.end()); // Erase duplicaes
        }


    }


};

class TEST_DBG : public ::testing::Test {
    protected:

    vector<string> reads;
    LL k;
    DBG_Reference_Implementation ref;
    BOSS<sdsl::bit_vector> boss;
    DBG dbg;

    void SetUp() override {
        reads = {"ATTCGTAGTCGTACGATAAATTTCGATGTAGGCTCGTTCGGTCGC", "ATTCTTTTCTTAGGCTAAAAAAAAA"};
        k = 3;
        ref = DBG_Reference_Implementation(reads,k);
        boss = build_BOSS_with_maps(reads, k, false);
        dbg = DBG(&boss);
    }

};


TEST_F(TEST_DBG, iterate_all_nodes){
    LL kmer_idx = 0;
    for(DBG::Node v : dbg.all_nodes()){
        string fetched = boss.get_node_label(v.id);
        string correct = ref.colex_kmers[kmer_idx++];
        cout << v.id << " " << fetched << endl;
        ASSERT_EQ(fetched, correct);
    }
    ASSERT_EQ(kmer_idx, ref.colex_kmers.size());
}

TEST_F(TEST_DBG, inedges){
    for(DBG::Node v : dbg.all_nodes()){
        LL in_idx = 0;
        for(DBG::Edge e : dbg.inedges(v)){
            ASSERT_EQ(e.dest, v.id);
            string kmer_from = boss.get_node_label(e.source);
            string kmer_to = boss.get_node_label(e.dest);
            ASSERT_EQ(kmer_to.back(), e.label);
            ASSERT_EQ(kmer_from, ref.inedges[kmer_to][in_idx++]);
        }

    }
}

TEST_F(TEST_DBG, outedges){
    for(DBG::Node v : dbg.all_nodes()){
        LL out_idx = 0;
        for(DBG::Edge e : dbg.outedges(v)){
            ASSERT_EQ(e.source, v.id);
            string kmer_from = boss.get_node_label(e.source);
            string kmer_to = boss.get_node_label(e.dest);
            ASSERT_EQ(kmer_to.back(), e.label);
            ASSERT_EQ(kmer_to, ref.outedges[kmer_from][out_idx++]);
        }
    }
}