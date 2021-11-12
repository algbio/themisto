#pragma once

#include <gtest/gtest.h>
#include <vector>
#include <unordered_map>
#include "setup_tests.hh"
#include "../globals.hh"
#include "../test_tools.hh"
#include "../libwheeler/BOSS.hh"
#include "../Themisto.hh"

// Testcase: put in a couple of reference sequences, sweep different k. For each k-mer, 
// ask what is the color set of that k-mer. It should coincide with the reads that contain
// that k-mer

struct ColoringTestCase{
    vector<string> references; //
    vector<string> colex_kmers; //
    unordered_map<string,set<LL> > kmer_to_ref_ids; //
    vector<set<LL> > color_sets; // kmer id -> color ids
    vector<LL> seq_id_to_color_id; //
    string fasta_data; //
    LL k; //
};

ColoringTestCase generate_testcase(vector<string> refs, vector<LL> colors, LL k){
    ColoringTestCase tcase;
    tcase.k = k;
    tcase.references = refs;
    set<string> kmer_set;
    for(LL i = 0; i < refs.size(); i++){
        tcase.fasta_data += ">\n" + refs[i] + "\n";
        tcase.seq_id_to_color_id.push_back(colors[i]);
        for(string kmer : get_all_distinct_kmers(refs[i], k)) kmer_set.insert(kmer);
    }
    vector<string> colex_kmers(kmer_set.begin(), kmer_set.end());
    sort(colex_kmers.begin(), colex_kmers.end(), colex_compare);
    tcase.colex_kmers = colex_kmers;
    vector<set<string> > ref_to_kmers;
    tcase.color_sets.resize(colex_kmers.size());

    // For all refs
    for(LL ref = 0; ref < tcase.references.size(); ref++){
        set<string> ref_kmers = get_all_distinct_kmers(tcase.references[ref], k);
        LL kmer_id = 0;

        // For all kmers of the whole data in colex order
        for(string kmer : tcase.colex_kmers){
            if(ref_kmers.count(kmer)){ // If this k-mer is found in the current reference, add the color
                tcase.color_sets[kmer_id].insert(colors[ref]);
                tcase.kmer_to_ref_ids[kmer].insert(ref);
            }
            kmer_id++;
        }
    }
    return tcase;
}


// All colors distinct
ColoringTestCase generate_testcase(vector<string> refs, LL k){
    vector<LL> colors;
    for(LL i = 0; i < refs.size(); i++){
        colors.push_back(i);
    }
    return generate_testcase(refs,colors,k);
}

vector<ColoringTestCase> generate_testcases(){
    vector<ColoringTestCase> cases;
    for(LL rep = 0; rep < 20; rep++){
        for(LL k = 1; k <= 20; k++){
            vector<string> refs;
            for(LL i = 0; i < 11; i++){
                if(rep % 2 == 0 && i == 10) // Add a duplicate to get redundant nodes
                    refs.push_back(refs.back());
                else
                    refs.push_back(get_random_dna_string(30, 2));
            }
            cases.push_back(generate_testcase(refs, k));
        }
    }
    logger << "Done generating testcases" << endl;
    return cases;
}

TEST(COLORING_TESTS, kmer_edge_case){
    // The coloring should not have a k-mer if it is not inside some (k+1)-mer
    vector<string> seqs = {"AACA", "ATAACAGACT"};
    vector<LL> colors = {0,1};
    LL k = 4;
    string fastafile = temp_file_manager.get_temp_file_name("");
    write_as_fasta(seqs, fastafile);
    BOSS<sdsl::bit_vector> boss = build_BOSS_with_maps(seqs, k, false);
    Coloring coloring;
    coloring.add_colors(boss, fastafile, colors, 1e6, 1, 1);
    
    LL node_id = boss.find_kmer("AACA");
    vector<LL> correct_colorset = {1}; // no color 0 because it does not have a (4+1)-mer
    set<LL> colorset_set = coloring.get_colorset(node_id, boss);
    vector<LL> colorset_vec(colorset_set.begin(), colorset_set.end());
    logger << "== Colorset edge case test == " << endl << correct_colorset << endl << colorset_vec << endl << "==" << endl;
    ASSERT_EQ(correct_colorset, colorset_vec);
}

TEST(COLORING_TESTS, random_testcases){
    for(ColoringTestCase tcase : generate_testcases()){
        logger << "Running testcase" << endl;
        string fastafilename = temp_file_manager.get_temp_file_name("ctest");
        throwing_ofstream fastafile(fastafilename);
        fastafile << tcase.fasta_data;
        fastafile.close();
        BOSS<sdsl::bit_vector> boss = build_BOSS_with_maps(tcase.references, tcase.k, false);
        Coloring coloring;
        coloring.add_colors(boss, fastafilename, tcase.seq_id_to_color_id, 1000, 3, 10);

        for(LL kmer_id = 0; kmer_id < tcase.colex_kmers.size(); kmer_id++){
            string kmer = tcase.colex_kmers[kmer_id];
            LL node_id = boss.find_kmer(kmer);
            set<LL> correct_colorset = tcase.color_sets[kmer_id];
            set<LL> colorset = coloring.get_colorset(node_id, boss);
            logger << colorset << endl << correct_colorset << endl;
            ASSERT_EQ(correct_colorset, colorset);
        }
    }
}
