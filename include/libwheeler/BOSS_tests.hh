#pragma once

#include "BOSS.hh"
#include "BOSS_builder.hh"
#include "../KMC_wrapper.hh"
#include "../tests/setup_tests.hh"
#include <gtest/gtest.h>

typedef BOSS<sdsl::bit_vector> boss_t;

class BOSS_TestCase{
    public:
    vector<string> reads;
    set<string> kmers;
    set<string> k_plus1_mers;
    set<char> alphabet;
    vector<string> colex_kmers;
    LL k;
};

// Test fixture
class BOSS_TEST : public ::testing::Test {
    public:
    vector<BOSS_TestCase> testcases;

    void SetUp() override {
        LL random_seed = 123674;
        srand(random_seed);
        LL n_reads = 10;
        for(float k_float = 1; k_float <= 255; k_float *= 1.5){ // kmer class supports at most k == 255
            LL k = k_float; // Cast to integer (floor)
            for(LL read_length = 1; read_length <= 128; read_length *= 2){
                k = min(k,(LL)254); // 254 + 2 = 256 is the biggest the EM sort can do
                BOSS_TestCase tcase;
                vector<string> reads;

                for(LL i = 0; i < n_reads; i++) 
                    reads.push_back(get_random_dna_string(read_length, 2 + rand()%2));

                tcase.reads = reads;

                for(string S : reads){
                    for(LL i = 0; i < (LL)S.size()-k; i++){
                        tcase.k_plus1_mers.insert(S.substr(i,k+1));
                        tcase.kmers.insert(S.substr(0,k));
                        tcase.kmers.insert(S.substr(1,k));
                    }
                }

                tcase.colex_kmers = vector<string>(tcase.kmers.begin(), tcase.kmers.end());
                sort(tcase.colex_kmers.begin(), tcase.colex_kmers.end(), colex_compare);

                for(string S : reads) for(char c : S) tcase.alphabet.insert(c);
                tcase.k = k;
                testcases.push_back(tcase);
            }
        }
    }
};

string BOSS_to_string(boss_t& boss){
    string bwt = boss.get_outlabels();
    string bwt_spaced;

    LL next = 0;
    for(LL i = 0; i < boss.outdegs_size(); i++){
        if(boss.outdegs_at(i) == 1) bwt_spaced += " ";
        else{
            bwt_spaced += bwt[next];
            next++;
        }
    }
    
    stringstream ss;
    ss << bwt_spaced << "\n" << boss.get_outdegs() << "\n" << boss.get_indegs() << "\n";
    ss << "C = ";
    vector<LL> C = boss.get_C_array();
    for(LL i = 0; i < (LL)C.size(); i++){
        if(i == 0 || C[i-1] != C[i]) ss << C[i] << " "; 
    }
    return ss.str();
}

void check_data_is_equal(boss_t& boss1, boss_t& boss2){
    ASSERT_TRUE(boss1 == boss2);
}

TEST_F(BOSS_TEST, serialization){
    
    for(BOSS_TestCase& tcase : testcases){
        boss_t boss = build_BOSS_with_maps(tcase.reads, tcase.k, false);

        string filename = get_temp_file_manager().create_filename();
        throwing_ofstream out(filename, ios::binary);
        boss.serialize(out.stream);
        out.stream.close();

        boss_t boss2;
        throwing_ifstream in(filename, ios::binary);
        boss2.load(in.stream);
        check_data_is_equal(boss, boss2);
    }
}

TEST_F(BOSS_TEST, search){
    for(BOSS_TestCase& tcase : testcases){
        // 1) If a k-mer is in set, it must be found
        // 2) If a k-mer is not in set, it must not be found
        boss_t boss = build_BOSS_with_maps(tcase.reads, tcase.k, false);

        // Check all the k-mers
        for(string kmer : tcase.kmers){
            LL boss_id_searched = boss.find_kmer(kmer);
            ASSERT_TRUE(boss_id_searched != -1);
        }

        // Try to check k-mers that are not found
        for(LL rep = 0; rep <= 1000; rep++){
            string kmer = get_random_string(tcase.k,20);
            if(tcase.kmers.count(kmer) == 0){
                // Not found
                ASSERT_TRUE(boss.find_kmer(kmer) == -1);
            }
        }
    }
}

TEST_F(BOSS_TEST, walk){
    for(BOSS_TestCase& tcase : testcases){
        boss_t boss = build_BOSS_with_maps(tcase.reads, tcase.k, false);
        ASSERT_TRUE(boss.walk(-1,*tcase.alphabet.begin()) == -1); // Walk from -1 gives -1
        for(string kmer : tcase.kmers){
            LL boss_id_searched = boss.find_kmer(kmer);
            ASSERT_TRUE(boss_id_searched != -1);
            for(char c : tcase.alphabet){
                bool found = tcase.k_plus1_mers.count(kmer + c);
                LL boss_id = boss.walk(boss_id_searched, c);
                if(found){
                    ASSERT_TRUE(boss_id != -1);
                    ASSERT_TRUE(kmer.substr(1) + c == boss.get_node_label(boss_id));
                }
                else ASSERT_TRUE(boss_id == -1);
            }
        }
    }
}

void test_construction(BOSS_TestCase& tcase, bool reverse_complements){
    boss_t boss_maps = build_BOSS_with_maps(tcase.reads, tcase.k, reverse_complements);

    // Build from KMC
    string fastafile = get_temp_file_manager().create_filename();
    throwing_ofstream out(fastafile);
    for(string S : tcase.reads) out << ">\n" << S << "\n";
    out.flush();
    string KMC_db_path_prefix = get_temp_file_manager().create_filename("KMC");
    KMC_wrapper(tcase.k+1, 1, 2, fastafile, get_temp_file_manager().get_dir(), KMC_db_path_prefix, reverse_complements, get_log_level() == LogLevel::OFF);
    Kmer_stream_from_KMC_DB kmer_stream(KMC_db_path_prefix, reverse_complements);
    BOSS_builder<boss_t, Kmer_stream_from_KMC_DB> bb;
    boss_t boss_KMC = bb.build(kmer_stream, 1e9, 2);

    set<string> A = boss_maps.get_all_edgemers();
    set<string> B = boss_KMC.get_all_edgemers();
    ASSERT_TRUE(A == B);
}


TEST_F(BOSS_TEST, construction){
    for(BOSS_TestCase& tcase : testcases){
        test_construction(tcase, false); // No reverse complements
        test_construction(tcase, true); // Reverse complements
    }
}

TEST_F(BOSS_TEST, node_label){
    for(BOSS_TestCase& tcase : testcases){
        boss_t boss = build_BOSS_with_maps(tcase.reads, tcase.k, false);
        for(string x : tcase.kmers){
            ASSERT_TRUE(boss.get_node_label(boss.find_kmer(x)) == x);
        }
    }
}


TEST_F(BOSS_TEST, empty){
    auto test = [](int k){
        boss_t boss(k);
        ASSERT_TRUE(boss.get_outlabels() == string(""));
        ASSERT_TRUE(boss.search_from_root("") == 0);
        ASSERT_TRUE(boss.search_from_root("A") == -1);
        ASSERT_TRUE(boss.search_from_root("") == 0);
        ASSERT_TRUE(boss.search_from_root("A") == -1);
        ASSERT_TRUE(boss.search_from_root("", 0) == 0);
        ASSERT_TRUE(boss.search_from_root("A", 1) == -1);
        ASSERT_TRUE(boss.walk(0, 'A') == -1);

        ASSERT_TRUE(boss.node_outlabels(0).size() == 0);

        ASSERT_TRUE(boss.indegree(0) == 0);
        ASSERT_TRUE(boss.outdegree(0) == 0);

        char label[256];
        boss.get_node_label(0, label);
        ASSERT_TRUE(label[0] == '\0');

        string kmer(k,'A');
        if(k == 0) ASSERT_TRUE(boss.find_kmer(kmer) == 0);
        else ASSERT_TRUE(boss.find_kmer(kmer) == -1);
    };

    test(1);
    test(10);
}