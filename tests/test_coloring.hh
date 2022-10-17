#pragma once

#include <gtest/gtest.h>
#include <vector>
#include <unordered_map>
#include "setup_tests.hh"
#include "globals.hh"
#include "test_tools.hh"
#include "sbwt/SBWT.hh"
#include "sbwt/globals.hh"
#include "new_coloring.hh"

// Testcase: put in a couple of reference sequences, sweep different k. For each k-mer,
// ask what is the color set of that k-mer. It should coincide with the reads that contain
// that k-mer

struct ColoringTestCase{
    vector<string> references; //
    vector<string> colex_kmers; //
    unordered_map<string,set<LL> > kmer_to_ref_ids; //
    vector<set<LL> > color_sets; // kmer id -> color ids
    vector<int64_t> seq_id_to_color_id; //
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

TEST(COLORING_TESTS, random_testcases){
    int64_t testcase_id = 0;
    for(ColoringTestCase tcase : generate_testcases()){
        logger << "Running testcase" << endl;
        string fastafilename = get_temp_file_manager().create_filename("ctest",".fna");
        sbwt::throwing_ofstream fastafile(fastafilename);
        fastafile << tcase.fasta_data;
        fastafile.close();
        plain_matrix_sbwt_t SBWT;
        build_nodeboss_in_memory<plain_matrix_sbwt_t>(tcase.references, SBWT, tcase.k, true);

        Coloring coloring;
        coloring.add_colors(SBWT, fastafilename, tcase.seq_id_to_color_id, 2048, 3, rand() % 3);

        for(LL kmer_id = 0; kmer_id < tcase.colex_kmers.size(); kmer_id++){
            string kmer = tcase.colex_kmers[kmer_id];
            LL node_id = SBWT.search(kmer);
            set<LL> correct_colorset = tcase.color_sets[kmer_id];
            vector<int64_t> colorvec = coloring.get_color_set_of_node_as_vector(node_id);
            set<LL> colorset(colorvec.begin(), colorvec.end());
            logger << node_id << ": " << colorset << " - " << correct_colorset << endl;
            ASSERT_EQ(correct_colorset, colorset);
        }
        testcase_id++;
    }
}

bool is_valid_kmer(const string& S){
    for(char c : S) if(c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
    return true;
}

template<typename color_set_t> requires Color_Set_Interface<color_set_t>
void test_coloring_on_coli3(plain_matrix_sbwt_t& matrix, string filename, std::vector<std::string>& seqs, int64_t k){

    std::vector<std::int64_t> colors;
    for(std::int64_t i = 0; i < seqs.size(); ++i){
        colors.push_back(i);
    }

    Coloring<color_set_t> c;
    c.add_colors(matrix, filename, colors, 1<<30, 3, 3);

    std::size_t seq_id = 0;
    write_log("Checking colors", LogLevel::MAJOR);
    for (const auto& seq : seqs) {

        std::cout << "Checking sequence: " << seq_id << "\n";
        //cout << seq << endl;

        const auto nodes = matrix.streaming_search(seq);

        for(LL i = 0; i < (LL)seq.size()-k+1; i++){
            LL node = nodes[i];
            if (node <= 0) {
                // This can happen if the input file has a non-ACGT character, as it does in this case
                ASSERT_FALSE(is_valid_kmer(seq.substr(i,k)));
            } else {
                auto vec = c.get_color_set_of_node_as_vector(node);
                const auto res = std::find(vec.begin(), vec.end(), seq_id);
                ASSERT_TRUE(res != vec.end());
            }
        }
        ++seq_id;
    }

    //TODO: also test that there are no extra colors in the color sets.
}

TEST(COLORING_TESTS, coli3) {
    std::string filename = "example_input/coli3.fna";

    std::vector<std::string> seqs;
    SeqIO::Unbuffered_Reader sr(filename);
    while(!sr.done()){
        string S = sr.get_next_query_stream().get_all();
        seqs.push_back(S);
    }

    const std::size_t k = 31;

    plain_matrix_sbwt_t::BuildConfig config;
    config.input_files = {filename};
    config.k = k;
    config.build_streaming_support = true;
    config.ram_gigas = 2;
    config.n_threads = 2;
    config.min_abundance = 1;
    plain_matrix_sbwt_t matrix(config);

    test_coloring_on_coli3<Bitmap_Or_Deltas_ColorSet>(matrix, filename, seqs, k);
    test_coloring_on_coli3<Fixed_Width_Int_Color_Set>(matrix, filename, seqs, k);
    test_coloring_on_coli3<Roaring_Color_Set>(matrix, filename, seqs, k);
}
