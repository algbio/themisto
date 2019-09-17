#pragma once

#include "BOSS.hh"
#include "TempFileManager.hh"
#include "BD_BWT_index.hh"
#include <algorithm>

class Coloring{

// The colors have names given by the user.
// Internally the colors are represented as integers.
// Each sequence corresponds to exactly one color, but one color
// can correspond to multiple sequence names

// sequence id = the rank of the read in the fasta file (starting from 0)
// sequence name = the fasta header of the sequence
// color name = the color names in the user-provided colorfile
// color id = the integer representation of a color name
// The color names are mapped to non-negative integers in lexicographic order

// The input for indexing is:
// - A fasta file with m sequences
// - A colorfile with m lines such that line i gives the
//   name of the color of the i-th sequence

// The input for querying is:
// - A fasta file of reads 
// The output for querying is:
// - For each read, a set of color names (not color ids)

public:

    //vector<set<color_t> > colors; // node id -> color set.

private:

    // private copy construtor and assignment operator: no copying, 
    // else might invalidate the pointer in the rank/select supports
    Coloring(const Coloring&); 
    Coloring& operator=(const Coloring&); 

public:

    sdsl::int_vector<> color_sets;
    sdsl::int_vector<> node_to_color_set_id;
    sdsl::bit_vector nonempty;
    sdsl::bit_vector color_set_starts;
    sdsl::rank_support_v<1> nonempty_rs;
    sdsl::select_support_mcl<1> color_set_starts_ss;
    sdsl::bit_vector nonempty_and_nonredundant;
    sdsl::bit_vector redundancy_marks; 
    sdsl::rank_support_v<1> nonempty_and_nonredundant_rs;
    LL n_colors;

    Coloring() {}

    // fastafile: the file containing the reference sequences in FASTA format
    // colors_assignments: for each sequence in the fastafile, an integer representing
    // the color of that sequence. The distinct colors should be a contiguous range of
    // integers from 0 to (number of colors -1).
    void add_colors(BOSS& boss, string fastafile, vector<LL> colors_assignments, LL k){
        n_colors = *std::max_element(colors_assignments.begin(), colors_assignments.end()) + 1;
        assert(is_contiguous_integer_range(colors_assignments, 0,n_colors-1));
        vector<pair<LL,LL> > node_color_pairs = get_node_color_pairs(boss, fastafile, colors_assignments, k);
        mark_redundant_color_sets(fastafile, boss);
        build_packed_representation(node_color_pairs, boss);
    }

    bool is_redundant(LL node){
        return redundancy_marks[node];
    }

    string to_string_internals(){
        stringstream ss;
        ss << "color_sets: " << color_sets << "\n";
        ss << "color_set_starts: " << color_set_starts << "\n";
        ss << "node_to_color_set_id: " << node_to_color_set_id << "\n";
        ss << "nonempty: " << nonempty << "\n";
        ss << "nonempty_and_nonredundant: " << nonempty_and_nonredundant << "\n";
        ss << "redundancy_marks: " << redundancy_marks << "\n";
        return ss.str();
    }

    string to_string(BOSS& boss){
        stringstream ss;
        for(LL node = 0; node < nonempty.size(); node++){
            ss << node << ": ";
            for(LL color : get_colorset(node, boss)) ss << color << " ";
            ss << "\n";
        }
        
        return ss.str();
    }

    // For debugging mostly
    vector<set<LL> > get_all_colorsets(BOSS& boss){
        vector<set<LL> > res;
        for(LL i = 0; i < nonempty.size(); i++){
            res.push_back(get_colorset(i, boss));
        }
        return res;
    }

    // Buffer must have enough space to accommodate all colors.
    // Puts the colors into the buffer starting from index 0 and
    // returns the number of colors added.
    // todo: tests
    LL get_colorset_to_buffer(LL node, BOSS& boss, vector<LL>& buffer){
        while(redundancy_marks[node] == 1) node = boss.get_predecessor(node); // Go to next non-redundant

        // The nodes that have an explicit color set are those that have a non-empty
        // colorset and are non-redundant.
        if(nonempty[node] == 0) {
            return 0;
        }
        else{
            LL rank = nonempty_and_nonredundant_rs.rank(node);
            LL color_set_id = node_to_color_set_id[rank];
            LL start = color_set_starts_ss.select(color_set_id+1);
            LL pos = start; // Position in internal color array
            LL idx = 0; // Position in outpuf buffer
            while(true){
                buffer[idx] = color_sets[pos];
                idx++;

                if(pos == color_sets.size()-1 || color_set_starts[pos+1] == 1) break;
                pos++;
            }
            return idx;
        }
    }

    // Super inefficient because allocates a new buffer every time
    set<LL> get_colorset(LL node, BOSS& boss){
        vector<LL> buffer(n_colors);
        LL count = get_colorset_to_buffer(node, boss, buffer);
        set<LL> result(buffer.begin(), buffer.begin() + count);
        return result;
    }

    void load_from_disk(string path_prefix){

        sdsl::load_from_file(color_sets, path_prefix + "colorsets");
        sdsl::load_from_file(node_to_color_set_id, path_prefix + "node-to-color-set-id");
        sdsl::load_from_file(nonempty, path_prefix + "colors-nonempty");
        sdsl::load_from_file(color_set_starts, path_prefix + "colors-starts");
        sdsl::load_from_file(redundancy_marks, path_prefix + "redundancy-marks");
        sdsl::load_from_file(nonempty_and_nonredundant, path_prefix + "nonempty-and-nonredundant");

        sdsl::load_from_file(nonempty_rs, path_prefix + "colors-nonempty-rs");
        nonempty_rs.set_vector(&nonempty);
        sdsl::load_from_file(color_set_starts_ss, path_prefix + "colors-starts-ss");
        color_set_starts_ss.set_vector(&color_set_starts);
        sdsl::load_from_file(nonempty_and_nonredundant_rs, path_prefix + "nonempty-and-nonredundant-rs");
        nonempty_and_nonredundant_rs.set_vector(&nonempty_and_nonredundant);

        n_colors = 0;
        for(LL i = 0; i < color_sets.size(); i++) n_colors = max(n_colors, (LL)(color_sets[i]+1));

    }

    void save_to_disk(string path_prefix){

        sdsl::store_to_file(color_sets, path_prefix + "colorsets");
        sdsl::store_to_file(node_to_color_set_id, path_prefix + "node-to-color-set-id");
        sdsl::store_to_file(nonempty, path_prefix + "colors-nonempty");
        sdsl::store_to_file(nonempty_rs, path_prefix + "colors-nonempty-rs");
        sdsl::store_to_file(color_set_starts, path_prefix + "colors-starts");
        sdsl::store_to_file(color_set_starts_ss, path_prefix + "colors-starts-ss");
        sdsl::store_to_file(redundancy_marks, path_prefix + "redundancy-marks");
        sdsl::store_to_file(nonempty_and_nonredundant, path_prefix + "nonempty-and-nonredundant");
        sdsl::store_to_file(nonempty_and_nonredundant_rs, path_prefix + "nonempty-and-nonredundant-rs");

    }

private:

    bool is_contiguous_integer_range(vector<LL>& v, LL from, LL to){
        set<LL> S(v.begin(), v.end());
        return *S.begin() >= from && *S.rbegin() <= to && S.size() == to - from + 1;
    }

    vector<pair<LL,LL> > get_node_color_pairs(BOSS& boss, string fastafile, vector<LL>& seq_id_to_color_id, LL k){
        write_log("Getting (node,color) pairs");
        vector<pair<LL,LL> > node_color_pairs;
        FASTA_reader fr(fastafile);
        LL n_processed = 0;
        while(!fr.done()){
            LL color = seq_id_to_color_id[n_processed];
            Read_stream rs = fr.get_next_query_stream();
            string seq = rs.get_all();
            write_log("Adding colors for sequence " + std::to_string(n_processed));
            if(seq.size() >= k){
                LL node = boss.search(seq.substr(0,k));
                node_color_pairs.push_back({node,color});
                //color_sets[node].insert(color);
                for(LL i = k; i < seq.size(); i++){
                    node = boss.walk(node, seq[i]);
                    assert(node != -1);
                    //color_sets[node].insert(color);
                    node_color_pairs.push_back({node,color});
                }
            }
            n_processed++;
        }

        node_color_pairs.shrink_to_fit();
        write_log("Sorting (node,color) pairs");
        std::sort(node_color_pairs.begin(), node_color_pairs.end());
        node_color_pairs.erase(std::unique(node_color_pairs.begin(), node_color_pairs.end()), 
                               node_color_pairs.end()); // Delete duplicates

        return node_color_pairs;
    }

    map<vector<LL>, vector<LL> > get_color_classes(const vector<pair<LL,LL> >& node_color_pairs){

        // keys: color classes, values: list of nodes with this color class
        map<vector<LL>, vector<LL> > class_to_nodes;

        LL run_start = 0;
        vector<LL> color_class;
        while(run_start < node_color_pairs.size()){
            LL run_end = run_start;
            LL node = node_color_pairs[run_start].first;
            while(run_end <= node_color_pairs.size()-2 && node_color_pairs[run_end+1].first == node){
                run_end++;
            }
            for(LL i = run_start; i <= run_end; i++){
                color_class.push_back(node_color_pairs[i].second);
            }
            class_to_nodes[color_class].push_back(node);
            color_class.clear();
            run_start = run_end + 1;
        }
        return class_to_nodes;
    }

    void build_packed_representation(const vector<pair<LL,LL> >& node_color_pairs, BOSS& boss){
        // todo: make these private functions return stuff instead of modifying the object state
        write_log("Building the packed representation");
        
        map<vector<LL>, vector<LL> > color_set_to_nodes = get_color_classes(node_color_pairs);

        LL n_classes = color_set_to_nodes.size();
        LL total_size_of_color_sets = 0;
        for(auto& P : color_set_to_nodes)
            total_size_of_color_sets += P.first.size();

        color_sets = sdsl::int_vector<>(total_size_of_color_sets, 0, ceil(log2(n_colors)));
        color_set_starts.resize(color_sets.size());
        sdsl::util::set_to_value(color_set_starts,0);

        nonempty.resize(boss.get_number_of_nodes());
        sdsl::util::set_to_value(nonempty,0);

        nonempty_and_nonredundant.resize(boss.get_number_of_nodes());
        sdsl::util::set_to_value(nonempty_and_nonredundant,0);

        node_to_color_set_id = sdsl::int_vector<>(boss.get_number_of_nodes(), 0, ceil(log2(n_classes)));

        LL index_in_color_sets = 0;
        LL class_id = 0;

        // Iterate all distinct color set
        for(auto& P : color_set_to_nodes){
            const vector<LL>& color_set = P.first;

            // Store the colors ids in the color set
            color_set_starts[index_in_color_sets] = 1;
            for(LL x : color_set){
                color_sets[index_in_color_sets] = x;
                index_in_color_sets++;
            }

            // Mark the nodes with this colorset to have a nonempty color set.
            for(LL node : P.second){
                nonempty[node] = 1;
                if(redundancy_marks[node] == 0){
                    nonempty_and_nonredundant[node] = 1;
                    node_to_color_set_id[node] = class_id;
                }
            }
            class_id++;
        }

        // Compactify node_to_color_set_id
        LL p = 0;
        for(LL i = 0; i < boss.get_number_of_nodes(); i++){
            if(nonempty_and_nonredundant[i]){
                node_to_color_set_id[p] = node_to_color_set_id[i];
                p++;
            }
        }
        node_to_color_set_id.resize(p);

        sdsl::util::init_support(nonempty_rs, &nonempty);
        sdsl::util::init_support(color_set_starts_ss, &color_set_starts);
        sdsl::util::init_support(nonempty_and_nonredundant_rs, &nonempty_and_nonredundant);
    }

    // Node color pairs should be sorted before calling this
    void mark_redundant_color_sets(string fastafile, BOSS& boss){

        // a node v is redundant if the color set of v must be necessarily
        // the same as the color set of pred(v). When can the color set change?
        //   (i) if some reference sequence ends at pred(v) (a color is lost)
        //   (ii) if some reference sequence starts at v (a color is gained)
        //   (iii) if pred(v) is branching (colors are split between siblings)
        //   (iv) if v has multiple predecessors (colors are gained from siblings)
        // We define that v is redundant if none of (i), (ii), (iii), (iv) hold.
        // If a node is redundant, then we do not need to intersect its colorset during
        // the pseudoalignment process. We can find the colorset of a redundant node by
        // going *backwards* through the unitig until we find a non-redundant node.

        // todo: this is expensive because it scans the whole reference data
        vector<string> first_and_last_kmers = get_first_and_last_kmers(fastafile, boss.get_k());

        write_log("Marking redundant color sets");

        redundancy_marks.resize(boss.get_number_of_nodes());
        sdsl::util::set_to_value(redundancy_marks,1); // 1 = redundant, 0 = non-redundant

        for(string kmer : first_and_last_kmers){
            LL id = boss.search(kmer);
            assert(id != -1);

            // cases (i) and (ii). Marking too much for quick testing.
            redundancy_marks[id] = 0; // case (ii): ref starts at id
            redundancy_marks[boss.get_successor(id)] = 0;  // case (i): ref ends at pred(id)

            // By the way this always marks at least one node, which is good beacuse if the
            // graph is just a single cycle (i.e. no branching nodes), then we need to mark
            // at least someting or else there could be no marks and then get_colorset() would
            // run forever
        }

        vector<char> buffer(256);
        LL mark_freq = (LL)log2(this->n_colors);
        for(LL node = 0; node < boss.get_number_of_nodes(); node++){
            if(boss.get_indegree(node) >= 2) redundancy_marks[node] = 0; // case (iv)
            else {
                LL pred = boss.get_predecessor(node);
                if(boss.is_branching(pred)) redundancy_marks[node] = 0; // case (iii)
            }

            if(boss.is_branching(node)){
                // Mark every log(n_colors)-th node as non-redundant on every unitig starting from here
                LL n_outedges = boss.get_outlabels_from(node, buffer);
                for(LL i = 0; i < n_outedges; i++){
                    LL next = boss.walk(node, buffer[i]);
                    LL distance = 0;
                    while(!boss.is_branching(next)){
                        // Will eventually find a branching node because we started from a branching
                        // node, so the graph can not be just a single cycle
                        next = boss.get_successor(next);
                        distance++;
                        if(distance == mark_freq){
                            redundancy_marks[next] = 0;
                            distance = 0;
                        }
                    }
                }
            }
        }

       //sdsl::util::set_to_value(redundancy_marks,0); // 1 = redundant, 0 = non-redundant. DEBUG!!!!!!

    }

};


// Testcase: put in a couple of reference sequences, sweep different k. For each k-mer, 
// ask what is the color set of that k-mer. It should coincide with the reads that contain
// that k-mer

class ColoringTester{
public:

    struct TestCase{
        vector<string> references; //
        vector<string> colex_kmers; //
        unordered_map<string,set<LL> > kmer_to_ref_ids; //
        vector<set<LL> > color_sets; // kmer id -> color ids
        vector<LL> seq_id_to_color_id; //
        string concat; //
        string fasta_data; //
        LL k; //
    };

    TestCase generate_testcase(vector<string> refs, vector<LL> colors, LL k){
        TestCase tcase;
        tcase.k = k;
        tcase.references = refs;
        for(LL i = 0; i < refs.size(); i++){    
            tcase.concat += tcase.references[i] + read_separator;
            tcase.seq_id_to_color_id.push_back(colors[i]);
            tcase.fasta_data += ">\n" + tcase.references[i] + "\n";
        }
        set<string> kmer_set = get_all_distinct_cyclic_kmers(tcase.concat + (char)(BD_BWT_index<>::END), k);
        vector<string> colex_kmers(kmer_set.begin(), kmer_set.end());
        sort(colex_kmers.begin(), colex_kmers.end(), colex_compare);
        tcase.colex_kmers = colex_kmers;
        vector<set<string> > ref_to_kmers;
        tcase.color_sets.resize(colex_kmers.size());
        for(LL ref = 0; ref < tcase.references.size(); ref++){
            set<string> ref_kmers = get_all_distinct_kmers(tcase.references[ref], k);
            LL kmer_id = 0;
            for(string kmer : tcase.colex_kmers){
                if(ref_kmers.count(kmer)){
                    tcase.color_sets[kmer_id].insert(colors[ref]);
                    tcase.kmer_to_ref_ids[kmer].insert(ref);
                }
                kmer_id++;
            }
        }
        return tcase;
    }

    // All colors distinct
    TestCase generate_testcase(vector<string> refs, LL k){
        vector<LL> colors;
        for(LL i = 0; i < refs.size(); i++){
            colors.push_back(i);
        }
        return generate_testcase(refs,colors,k);
    }

    vector<TestCase> generate_testcases(){
        vector<TestCase> cases;
        for(LL rep = 0; rep < 20; rep++){
            for(LL k = 1; k <= 20; k++){
                vector<string> refs;
                for(LL i = 0; i < 11; i++){
                    if(rep % 2 == 0 && i == 10) // Add a duplicate to get redundant nodes
                        refs.push_back(refs.back());
                    else
                        refs.push_back(get_random_string(30, 2));
                }
                cases.push_back(generate_testcase(refs, k));
            }
        }
        cerr << "Done generating testcases" << endl;
        return cases;
    }

    void run_testcase(TestCase tcase){
        string fastafilename = temp_file_manager.get_temp_file_name("ctest");
        ofstream fastafile(fastafilename);
        fastafile << tcase.fasta_data;
        fastafile.close();
        BOSS boss = build_BOSS_with_bibwt(tcase.concat, tcase.k);
        Coloring coloring;
        coloring.add_colors(boss, fastafilename, tcase.seq_id_to_color_id, tcase.k);

        cout << tcase.colex_kmers << endl;
        cout << coloring.to_string_internals() << endl;

        /*cout << tcase.references << endl;
        cout << coloring.to_string(boss) << endl;
        cout << boss.to_string() << endl;
        for(LL i = 0; i < tcase.colex_kmers.size(); i++){
            cout << i << " " << tcase.colex_kmers[i] << endl;
        }
        cout << "True color sets: " << endl;
        for(LL i = 0; i < tcase.colex_kmers.size(); i++){
            cout << i << " " << tcase.colex_kmers[i] << " " << tcase.kmer_to_ref_ids[tcase.colex_kmers[i]] << endl;
        }*/
        for(LL node_id = 0; node_id < tcase.colex_kmers.size(); node_id++){
            //string kmer = tcase.colex_kmers[node_id];
            set<LL> correct_colorset = tcase.color_sets[node_id];
            set<LL> colorset = coloring.get_colorset(node_id, boss);
            assert(correct_colorset == colorset);
        }
    }

    void test(){
        for(TestCase tcase : generate_testcases()){
            run_testcase(tcase);
        }
    }

};

void test_coloring(){
    ColoringTester tester;
    tester.test();
}