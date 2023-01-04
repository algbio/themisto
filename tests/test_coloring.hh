#pragma once

#include <gtest/gtest.h>
#include <vector>
#include <unordered_map>
#include "setup_tests.hh"
#include "globals.hh"
#include "test_tools.hh"
#include "sbwt/SBWT.hh"
#include "sbwt/globals.hh"
#include "sbwt/throwing_streams.hh"
#include "extract_unitigs.hh"
#include "DBG.hh"
#include "coloring/Coloring.hh"
#include "coloring/Coloring_Builder.hh"

// Testcase: put in a couple of reference sequences, sweep different k. For each k-mer,
// ask what is the color set of that k-mer. It should coincide with the reads that contain
// that k-mer

struct ColoringTestCase{
    vector<string> references; //
    vector<string> colex_kmers; //
    unordered_map<string,set<int64_t> > kmer_to_ref_ids; //
    vector<set<int64_t> > color_sets; // kmer id -> color ids
    vector<int64_t> seq_id_to_color_id; //
    string fasta_data; //
    int64_t k; //
};

ColoringTestCase generate_testcase(vector<string> refs, vector<int64_t> colors, int64_t k){
    ColoringTestCase tcase;
    tcase.k = k;
    tcase.references = refs;
    set<string> kmer_set;
    for(int64_t i = 0; i < refs.size(); i++){
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
    for(int64_t ref = 0; ref < tcase.references.size(); ref++){
        set<string> ref_kmers = get_all_distinct_kmers(tcase.references[ref], k);
        int64_t kmer_id = 0;

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
ColoringTestCase generate_testcase(vector<string> refs, int64_t k){
    vector<int64_t> colors;
    for(int64_t i = 0; i < refs.size(); i++){
        colors.push_back(i);
    }
    return generate_testcase(refs,colors,k);
}

vector<ColoringTestCase> generate_testcases(){
    vector<ColoringTestCase> cases;
    for(int64_t rep = 0; rep < 20; rep++){
        for(int64_t k = 1; k <= 20; k++){
            vector<string> refs;
            for(int64_t i = 0; i < 11; i++){
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

        Coloring<> coloring;
        Coloring_Builder<> cb;
        sbwt::SeqIO::Reader<> reader(fastafilename);
        cb.build_coloring(coloring, SBWT, reader, tcase.seq_id_to_color_id, 2048, 3, rand() % 3);

        for(int64_t kmer_id = 0; kmer_id < tcase.colex_kmers.size(); kmer_id++){
            string kmer = tcase.colex_kmers[kmer_id];
            int64_t node_id = SBWT.search(kmer);
            set<int64_t> correct_colorset = tcase.color_sets[kmer_id];
            vector<int64_t> colorvec = coloring.get_color_set_of_node_as_vector(node_id);
            set<int64_t> colorset(colorvec.begin(), colorvec.end());
            logger << node_id << ": " << colorset << " - " << correct_colorset << endl;
            ASSERT_EQ(correct_colorset, colorset);
        }
        testcase_id++;
    }
}

bool is_valid_kmer(const char* S, int64_t k){
    for(int64_t i = 0; i < k; i++){
        char c = S[i];
        if(c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
    }
    return true;
}


bool is_valid_kmer(const string& S){
    for(char c : S) if(c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
    return true;
}

// filename should not have added reverse complements
// seqs and seq_to_color chould have added reverse complements
template<typename color_set_t, typename color_set_view_t> requires Color_Set_Interface<color_set_t>
void test_coloring_on_coli3(plain_matrix_sbwt_t& matrix, string filename, std::vector<std::string>& seqs, vector<int64_t>& seq_to_color, int64_t k){

    std::unordered_map<Kmer<32>, vector<int64_t>> true_colors; // True color sets of k-mers

    write_log("Hashing true color sets", LogLevel::MAJOR);
    for(int64_t i = 0; i < seqs.size(); ++i){
        for(int64_t j = 0; j < (int64_t)seqs[i].size() - k + 1; j++){
            if(is_valid_kmer(seqs[i].c_str() + j, k)){
                Kmer<32> x(seqs[i].c_str() + j, k);
                true_colors[x].push_back(seq_to_color[i]);
            }
        }
    }

    // Deduplicate colors in color sets
    for(auto& [kmer, colorset] : true_colors){
        std::sort(colorset.begin(), colorset.end());
        colorset.erase( std::unique( colorset.begin(), colorset.end() ), colorset.end() );
    }

    Coloring<color_set_t> c;
    Coloring_Builder<color_set_t> cb;
    sbwt::SeqIO::Reader reader(filename); reader.enable_reverse_complements();
    cb.build_coloring(c, matrix, reader, seq_to_color, 1<<30, 3, 3);

    write_log("Checking colors", LogLevel::MAJOR);
    DBG dbg(&matrix);
    for(DBG::Node v : dbg.all_nodes()){
        
        vector<int64_t> colorset = c.get_color_set_of_node_as_vector(v.id);
        std::sort(colorset.begin(), colorset.end());

        Kmer<32> kmer(dbg.get_node_label(v));
        const vector<int64_t>& ref = true_colors[kmer];

        //logger << kmer.to_string() << " " << vec_to_string(colorset) << " " << vec_to_string(ref) << endl;

        ASSERT_EQ(colorset, ref);
        
    }
}

// Seqs need to contain reverse complements
// Colors should contain added reverse complements
// The sequences in the fasta filename should *not* contain added reverse complements (it's for ggcat)
void test_construction_from_colored_unitigs(plain_matrix_sbwt_t& SBWT, const vector<string>& seqs, vector<int64_t> seq_to_color, string filename, int64_t k){

    Coloring<SDSL_Variant_Color_Set> coloring;
    Coloring_Builder<SDSL_Variant_Color_Set> cb;
    sbwt::SeqIO::Reader reader(filename); reader.enable_reverse_complements();
    cb.build_coloring(coloring, SBWT, reader, seq_to_color, 1<<30, 3, 3);

    UnitigExtractor<Coloring<SDSL_Variant_Color_Set>> UE;
    DBG dbg(&SBWT);

    string unitigs_outfile = get_temp_file_manager().create_filename("unitigs-",".fna");
    string unitig_colors_outfile = get_temp_file_manager().create_filename("unitigs-colors-",".txt");

    sbwt::throwing_ofstream unitigs_out(unitigs_outfile);
    sbwt::throwing_ofstream unitig_colors_out(unitig_colors_outfile);

    sbwt::SeqIO::NullStream gfa_null_stream;
    UE.extract_unitigs(dbg, coloring, unitigs_out.stream, true, unitig_colors_out.stream, gfa_null_stream, 0);

    unitigs_out.close();
    unitig_colors_out.close();

    // Parse unitigs and colors from disk
    vector<string> unitigs;
    vector<vector<int64_t> > color_sets;
    SeqIO::Reader<> unitigs_in(unitigs_outfile);
    while(true){ // Read unitigs
        string S = unitigs_in.get_next_read();
        if(S == "") break;
        else unitigs.push_back(S);
    }

    // Parse unitig colors
    sbwt::throwing_ifstream color_sets_in(unitig_colors_outfile);
    string line;
    while(getline(color_sets_in.stream, line)){ // Read colors
        vector<string> tokens = split(line);
        vector<int64_t> colors;
        for(int64_t i = 1; i < tokens.size(); i++){ // i == 0 is the unitig id, which we ignore here
            colors.push_back(fast_string_to_int(tokens[i].c_str(), tokens[i].size()));
        }
        color_sets.push_back(colors);
    }


    // Build from GGCAT
    // Split the unitigs into one unitig per file for ggcat
    vector<string> ggcat_input_files = split_seqs_to_separate_files(filename);
    Coloring<SDSL_Variant_Color_Set> coloring2;
    Coloring_Builder<SDSL_Variant_Color_Set> cb2;
    sbwt::SeqIO::Reader reader2(filename); reader2.enable_reverse_complements();
    GGCAT_unitig_database db(ggcat_input_files, 1, k, 3, true);
    cb2.build_from_colored_unitigs(coloring2, reader2, SBWT, 1<<30, 3, 3, db);

    // Compare

    ASSERT_EQ(coloring.largest_color(), coloring2.largest_color());

    for(DBG::Node v : dbg.all_nodes()){
        vector<int64_t> c1 = coloring.get_color_set_of_node(v.id).get_colors_as_vector();
        vector<int64_t> c2 = coloring2.get_color_set_of_node(v.id).get_colors_as_vector();
        ASSERT_EQ(c1, c2);
    }

}

TEST(COLORING_TESTS, coli3) {
    std::string filename = "example_input/coli3.fna";
    int64_t k = 31;

    // Build the SBWT
    string fastafile = get_temp_file_manager().create_filename("",".fna");
    string indexprefix = get_temp_file_manager().create_filename();
    string tempdir = get_temp_file_manager().get_dir();
    vector<string> args = {"build", "-k", to_string(k), "-i", filename, "-o", indexprefix, "--temp-dir", tempdir, "--reverse-complements"};
    cout << args << endl;
    sbwt::Argv argv(args);
    build_index_main(argv.size, argv.array);
    plain_matrix_sbwt_t SBWT;
    SBWT.load(indexprefix + ".tdbg");

    // Read the sequences
    SeqIO::Unbuffered_Reader sr(filename);
    vector<string> seqs;
    vector<int64_t> seq_to_color;
    int64_t seq_idx = 0;
    while(!sr.done()){
        string S = sr.get_next_query_stream().get_all();
        seqs.push_back(S);
        seqs.push_back(get_rc(S));
        seq_to_color.push_back(seq_idx);
        seq_to_color.push_back(seq_idx); // For reverse complements. Yeah this is getting real hacky.
        seq_idx++;
    }

    write_log("Testing construction from colored unitigs", LogLevel::MAJOR);
    test_construction_from_colored_unitigs(SBWT, seqs, seq_to_color, filename, k);

    write_log("Testing Standard color set", LogLevel::MAJOR);
    test_coloring_on_coli3<SDSL_Variant_Color_Set, SDSL_Variant_Color_Set_View>(SBWT, filename, seqs, seq_to_color, k);
    write_log("Testing Roaring_Color_Set", LogLevel::MAJOR);
    test_coloring_on_coli3<Roaring_Color_Set, Roaring_Color_Set>(SBWT, filename, seqs, seq_to_color, k);

}
