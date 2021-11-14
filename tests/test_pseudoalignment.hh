#pragma once

#include <iostream>
#include <vector>
#include <unordered_map>
#include "../globals.hh"
#include "../test_tools.hh"
#include "../Themisto.hh"
#include "setup_tests.hh"
#include <gtest/gtest.h>

using namespace std;

LL random_seed = 123674;

class TestCase{
    public:

    vector<string> genomes;
    unordered_map<string, set<LL> > node_to_color_ids; // kmer- > set of color ids
    vector<string> queries;
    vector<string> colex_kmers;
    LL n_colors; // Distinct colors
    LL k;

    vector<LL> seq_to_color_id; // for each sequence: the color of that sequence
};


vector<TestCase> generate_testcases(LL genome_length, LL n_genomes, LL n_queries, LL query_length, LL min_k, LL max_k, LL n_colors){
    srand(random_seed);
    vector<TestCase> testcases;

    for(LL rep = 0; rep < 5; rep++){
        for(LL k = min_k; k <= max_k; k++){
            TestCase tcase;
            tcase.k = k;

            // Build genomes and assign color ids
            for(LL i = 0; i < n_genomes; i++){
                tcase.genomes.push_back(get_random_dna_string(genome_length,2));
                tcase.seq_to_color_id.push_back(rand() % n_colors);
            }

            tcase.n_colors = *std::max_element(tcase.seq_to_color_id.begin(), tcase.seq_to_color_id.end()) + 1;

            // Get all k-mers and colex-sort them
            set<string> all_kmers;
            for(LL i = 0; i < n_genomes; i++){
                for(string kmer : get_all_distinct_kmers(tcase.genomes[i], k)) all_kmers.insert(kmer);
            }
            
            tcase.colex_kmers = vector<string>(all_kmers.begin(), all_kmers.end());
            sort(tcase.colex_kmers.begin(), tcase.colex_kmers.end(), colex_compare);

            // List k-mer sets for each color
            vector<set<string> > color_to_kmer_set(n_colors);
            for(LL genome_id = 0; genome_id < tcase.genomes.size(); genome_id++){
                string genome = tcase.genomes[genome_id];
                LL color_id = tcase.seq_to_color_id[genome_id];
                for(string kmer : get_all_distinct_kmers(genome,k)){
                    color_to_kmer_set[color_id].insert(kmer);
                }
            }

            // List all color names for each k-mer (!= each node)
            for(string kmer : tcase.colex_kmers){
                set<LL> colorset;
                for(LL color_id = 0; color_id < n_colors; color_id++){
                    if(color_to_kmer_set[color_id].count(kmer)){
                        colorset.insert(color_id);
                    }
                }
                tcase.node_to_color_ids[kmer] = colorset;
            }
            

            // Build queries
            for(LL i = 0; i < n_queries; i++) tcase.queries.push_back(get_random_dna_string(query_length,2));

            testcases.push_back(tcase);
        }
    }
    return testcases;
}

// Returns set of color names
set<LL> pseudoalign_to_colors_trivial(string& query, TestCase& tcase, bool reverse_complements){
    set<LL> alignments;
    for(LL i = 0; i < tcase.n_colors; i++) alignments.insert(i); // All color names

    bool at_least_one = false;
    // For each k-mer in query, get the color name set and intersect that with alignments
    for(string kmer : get_all_kmers(query, tcase.k)){
        set<LL> colorset = tcase.node_to_color_ids[kmer];
        if(reverse_complements){
            set<LL> colorset_rc = tcase.node_to_color_ids[get_rc(kmer)];
            for(LL color : colorset_rc) colorset.insert(color);
        }
        if(colorset.size() >= 1) {
            at_least_one = true;
            alignments = intersect(alignments, colorset);
        }
    }

    if(at_least_one == false) alignments.clear();
    return alignments;
}


TEST(TEST_PSEUDOALIGN, random_testcases){
    logger << "Testing pseudolign" << endl;

    LL testcase_id = -1;
    LL ref_length = 100;
    LL n_refs = 50;
    LL n_threads = 3;
    LL n_queries = 10000;
    LL query_length = 20;
    LL k_min = 1;
    LL k_max = 20;
    LL n_colors = 5;
    for(TestCase tcase : generate_testcases(ref_length, n_refs, n_queries, query_length, k_min,k_max, n_colors)){
        testcase_id++;
        logger << "Running alignment testcase" << endl;

        string genomes_outfilename = create_temp_filename("genomes-");
        string colors_outfilename = create_temp_filename("colors-");
        string queries_outfilename = create_temp_filename("queries-");
        string colorfile_outfilename = create_temp_filename("colorfile-");
        string boss_outfilename = create_temp_filename("boss-");
        string coloring_outfilename = create_temp_filename("coloring-");

        throwing_ofstream genomes_out(genomes_outfilename);
        for(string genome : tcase.genomes){
            genomes_out << ">\n" << genome << "\n";
        }
        genomes_out.close();

        throwing_ofstream colors_out(colorfile_outfilename);
        for(LL i = 0; i < tcase.seq_to_color_id.size(); i++){
            colors_out << tcase.seq_to_color_id[i] << "\n";
        }
        colors_out.close();

        throwing_ofstream queries_out(queries_outfilename);
        for(string query : tcase.queries){
            queries_out << ">\n" << query << "\n";
        }
        queries_out.close();

        {
            Themisto kl_build;
            kl_build.construct_boss(genomes_outfilename, tcase.k, 1000, 2, true);
            kl_build.construct_colors(genomes_outfilename, colorfile_outfilename, 1000, 3, 10);
            kl_build.save_boss(boss_outfilename);
            kl_build.save_colors(coloring_outfilename);
        }

        Themisto kl; // Load back to test serialization
        kl.load_boss(boss_outfilename);
        kl.load_colors(coloring_outfilename);

        // Check against reference boss
        BOSS<sdsl::bit_vector> ref_boss = build_BOSS_with_maps(tcase.genomes, tcase.k, true); // Reverse complements enabled because the KMC construction uses those
        ASSERT_TRUE(ref_boss == kl.boss);

        // Run without rc
        string final_file = create_temp_filename("finalfile-");
        Sequence_Reader sr(queries_outfilename, FASTA_MODE);
        kl.pseudoalign_parallel(n_threads, sr, final_file, false, 300, false, true);
        vector<set<LL> > our_results = kl.parse_output_format_from_disk(final_file);

        // Run with rc
        string final_file_rc = create_temp_filename("finalfile_rc-");
        Sequence_Reader sr2(queries_outfilename, FASTA_MODE);
        kl.pseudoalign_parallel(n_threads, sr2, final_file_rc, true, 300, false, true);
        vector<set<LL> > our_results_rc = kl.parse_output_format_from_disk(final_file_rc);

        for(LL i = 0; i < tcase.queries.size(); i++){
            string query = tcase.queries[i];

            set<LL> brute = pseudoalign_to_colors_trivial(query, tcase, false);
            set<LL> brute_rc = pseudoalign_to_colors_trivial(query, tcase, true);

            logger << brute << endl << brute_rc << "-" << endl;

            ASSERT_EQ(brute, our_results[i]);
            ASSERT_EQ(brute_rc, our_results_rc[i]);
        }
    }
}
