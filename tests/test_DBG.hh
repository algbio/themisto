#pragma once

#include "setup_tests.hh"
#include <gtest/gtest.h>
#include "DBG.hh"
#include "BOSS_builder.hh"

TEST(TEST_DBG, basic){
    vector<string> reads = {"ATTCGTAGTCGTACGAT", "ATTCTTTTCTTAGGCTAAAAAAAAA"};
    LL k = 4;

    vector<string> kmers;
    for(string S : reads) 
        for(string x : get_all_distinct_kmers(S, k)) 
            kmers.push_back(x);
    std::sort(kmers.begin(), kmers.end(), colex_compare);
    kmers.erase(unique(kmers.begin(), kmers.end()), kmers.end()); // Erase duplicaes
    cout << kmers << endl;
    
    BOSS<sdsl::bit_vector> boss = build_BOSS_with_maps(reads, k, false);
    DBG dbg(&boss);
    cout << dbg.is_dummy << endl;
    
    LL kmer_idx = 0;
    for(DBG::Node v : dbg.all_nodes()){
        string fetched = boss.get_node_label(v.id);
        string correct = kmers[kmer_idx++];
        cout << v.id << " " << fetched << " " << correct << endl;
        ASSERT_EQ(fetched, correct);
    }
}
