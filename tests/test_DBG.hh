#pragma once

#include "setup_tests.hh"
#include <gtest/gtest.h>
#include "DBG.hh"
#include "BOSS_builder.hh"

class DBG_Reference_Implementation{

public:

    vector<string> colex_kmers;
    unordered_map<string, vector<string>> outedges; // Outedges from a node are colex-sorted by destination

    DBG_Reference_Implementation(vector<string> reads, LL k){
        // Get sorted k-mers
        for(string S : reads) 
            if(S.size() >= k+1) // Edge centric
                for(string x : get_all_distinct_kmers(S, k)) 
                    colex_kmers.push_back(x);
        std::sort(colex_kmers.begin(), colex_kmers.end(), colex_compare);
        colex_kmers.erase(unique(colex_kmers.begin(), colex_kmers.end()), colex_kmers.end()); // Erase duplicaes
        cout << colex_kmers << endl;

        // Get outedges
        for(string S : reads){
            for(string x : get_all_distinct_kmers(S, k+1)){
                string source = x.substr(0,k);
                string dest = x.substr(1,k);
                outedges[source].push_back(dest);
            }
        }

        // Sort out-edges and delete duplicates
        for(auto& keyval : outedges){
            vector<string>& out = keyval.second;
            std::sort(out.begin(), out.end(), colex_compare);
            out.erase(unique(out.begin(), out.end()), out.end()); // Erase duplicaes
        }

    }


};

TEST(TEST_DBG, basic){
    vector<string> reads = {"ATTCGTAGTCGTACGATAAATTTCGATGTAGGCTCGTTCGGTCGC", "ATTCTTTTCTTAGGCTAAAAAAAAA"};
    LL k = 3;

    DBG_Reference_Implementation DBG_ref(reads, k);
    
    BOSS<sdsl::bit_vector> boss = build_BOSS_with_maps(reads, k, false);
    DBG dbg(&boss);
    cout << dbg.is_dummy << endl;
    
    LL kmer_idx = 0;
    for(DBG::Node v : dbg.all_nodes()){
        string fetched = boss.get_node_label(v.id);
        string correct = DBG_ref.colex_kmers[kmer_idx++];
        cout << v.id << " " << fetched << " " << correct << endl;
        ASSERT_EQ(fetched, correct);
        LL edge_idx = 0;
        for(DBG::Edge e : dbg.outedges(v)){
            cout << e.source << " -> " << e.dest << " " << e.label << endl;
            ASSERT_EQ(e.source, v.id);
            string kmer_from = boss.get_node_label(e.source);
            string kmer_to = boss.get_node_label(e.dest);
            ASSERT_EQ(kmer_to.back(), e.label);
            ASSERT_EQ(kmer_to, DBG_ref.outedges[kmer_from][edge_idx]);
            edge_idx++;
        }
        
    }
    ASSERT_EQ(kmer_idx, DBG_ref.colex_kmers.size());
}
