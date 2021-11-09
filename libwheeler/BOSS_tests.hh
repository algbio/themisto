#pragma once
#include "BOSS.hh"
#include "BOSS_builder.hh"

template<typename boss_t = BOSS<sdsl::bit_vector>>
class Boss_Tester{
public:

    LL random_seed = 123674;

    class TestCase{
        public:
        vector<string> reads;
        set<string> kmers;
        set<string> k_plus1_mers;
        set<char> alphabet;
        vector<string> colex_kmers;
        LL k;
    };

    Boss_Tester(){
        
    }

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
        assert(boss1 == boss2);
    }

    void test_serialization(TestCase& tcase){
        boss_t boss = build_BOSS_with_maps(tcase.reads, tcase.k, false);
        boss.save_to_disk("test_out/");
        boss_t boss2;
        boss2.load_from_disk("test_out/");
        check_data_is_equal(boss, boss2);
    }

    vector<LL> get_neighbors(TestCase& tcase, LL id){
        string kmer = tcase.colex_kmers[id];
        vector<LL> neighbors;
        for(char c : tcase.alphabet){
            if(tcase.k_plus1_mers.count(kmer + c)){
                neighbors.push_back(tcase.kmer_to_node_id[kmer.substr(1)+c]);
            }
        }
        return neighbors;
    }

    vector<Boss_Tester::TestCase> generate_testcases(){
        srand(random_seed);
        vector<TestCase> testcases;
        LL n_reads = 10;
        for(float k_float = 1; k_float <= 255; k_float *= 1.5){ // kmer class supports at most k == 255
            LL k = k_float; // Cast to integer (floor)
            for(LL read_length = 1; read_length <= 128; read_length *= 2){
                k = min(k,(LL)254); // 254 + 2 = 256 is the biggest the EM sort can do
                TestCase tcase;
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
        return testcases;
    }

    void test_search(TestCase& tcase){
        // 1) If a k-mer is in set, it must be found
        // 2) If a k-mer is not in set, it must not be found
        boss_t boss = build_BOSS_with_maps(tcase.reads, tcase.k, false);

        // Check all the k-mers
        for(string kmer : tcase.kmers){
            LL boss_id_searched = boss.find_kmer(kmer);
            assert(boss_id_searched != -1);
        }

        // Try to check k-mers that are not found
        for(LL rep = 0; rep <= 1000; rep++){
            string kmer = get_random_string(tcase.k,20);
            if(tcase.kmers.count(kmer) == 0){
                // Not found
                assert(boss.find_kmer(kmer) == -1);
            }
        }
    }

    void test_walk(TestCase& tcase){
        boss_t boss = build_BOSS_with_maps(tcase.reads, tcase.k, false);
        cout << boss.to_string() << endl;
        assert(boss.walk(-1,*tcase.alphabet.begin()) == -1); // Walk from -1 gives -1
        for(string kmer : tcase.kmers){
            cout << "Testing searching with " << kmer << endl;
            LL boss_id_searched = boss.find_kmer(kmer);
            assert(boss_id_searched != -1);
            for(char c : tcase.alphabet){
                bool found = tcase.k_plus1_mers.count(kmer + c);
                LL boss_id = boss.walk(boss_id_searched, c);
                if(found){
                    assert(boss_id != -1);
                    assert(kmer.substr(1) + c == boss.get_node_label(boss_id));
                }
                else assert(boss_id == -1);
            }
        }
    }

    void test_construction(TestCase& tcase, bool reverse_complements){
        boss_t boss_maps = build_BOSS_with_maps(tcase.reads, tcase.k, reverse_complements);

        // Build from KMC
        
        string fastafile = temp_file_manager.get_temp_file_name("");
        throwing_ofstream out(fastafile);
        for(string S : tcase.reads) out << ">\n" << S << "\n";
        out.flush();
        cout << "Calling KMC" << endl;
        string KMC_db = KMC_wrapper(tcase.k+1, 1, 2, fastafile, temp_file_manager.get_dir(), reverse_complements);
        cout << "KMC DONE" << endl;
        Kmer_stream_from_KMC_DB kmer_stream(KMC_db, reverse_complements);
        BOSS_builder<boss_t, Kmer_stream_from_KMC_DB> bb;
        boss_t boss_KMC = bb.build(kmer_stream, 1e9, 2);

        cout << "maps" << endl << boss_maps.to_string() << endl;
        cout << "kmc" << endl << boss_KMC.to_string() << endl;
        set<string> A = boss_maps.get_all_edgemers();
        set<string> B = boss_KMC.get_all_edgemers();
        assert(A == B);
        

    }

    void test_get_node_label(TestCase& tcase){
        boss_t boss = build_BOSS_with_maps(tcase.reads, tcase.k, false);
        for(string x : tcase.kmers){
            assert(boss.get_node_label(boss.find_kmer(x)) == x);
        }
    }
};

template<typename boss_t = BOSS<sdsl::bit_vector>>
void test_empty_BOSS(){

    auto test = [](int k){
        boss_t boss(k);
        assert(boss.get_outlabels() == string(""));
        assert(boss.search_from_root("") == 0);
        assert(boss.search_from_root("A") == -1);
        assert(boss.search_from_root("") == 0);
        assert(boss.search_from_root("A") == -1);
        assert(boss.search_from_root("", 0) == 0);
        assert(boss.search_from_root("A", 1) == -1);
        assert(boss.walk(0, 'A') == -1);

        assert(boss.node_outlabels(0).size() == 0);

        assert(boss.indegree(0) == 0);
        assert(boss.outdegree(0) == 0);

        char label[256];
        boss.get_node_label(0, label);
        assert(label[0] == '\0');

        string kmer(k,'A');
        if(k == 0) assert(boss.find_kmer(kmer) == 0);
        else assert(boss.find_kmer(kmer) == -1);
    };

    test(1);
    test(10);
}

template<typename boss_t = BOSS<sdsl::bit_vector>>
void test_BOSS(){
    cerr << "Testing BOSS" << endl;
    disable_logging();

    test_empty_BOSS();
    cerr << "Empty BOSS test passed" << endl;

    // Test self-assignment
    boss_t boss;
    boss = boss; // Should not crash

    Boss_Tester<boss_t> tester;
    for(typename Boss_Tester<boss_t>::TestCase tcase : tester.generate_testcases()){
        cerr << "Running testcase: " << "k = " << tcase.k << ", read_length = " << tcase.reads[0].size() << " n_reads = " << tcase.reads.size() << endl;
        cerr << "Testing construction with reverse complements" << endl;
        tester.test_construction(tcase, true);
        cerr << "Testing construction without reverse complements" << endl;
        tester.test_construction(tcase, false);
        cerr << "Testing serialization" << endl;
        tester.test_serialization(tcase);
        cerr << "Testing search" << endl;
        tester.test_search(tcase);
        cerr << "Testing walk" << endl;
        tester.test_walk(tcase);
        cerr << "Testing get_node_label" << endl;
        tester.test_get_node_label(tcase);
    }

    cerr << "BOSS tests passed" << endl;
    enable_logging();
}
