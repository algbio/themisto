#pragma once

#include "BOSS.hh"
#include "TempFileManager.hh"
#include "BD_BWT_index.hh"
#include "EM_algorithms.hh"
#include "WorkDispatcher.hh"
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

    class AlignerThread : public DispatcherConsumerCallback{

        private:

        AlignerThread(const AlignerThread&); // No copying
        AlignerThread& operator=(const AlignerThread&); // No copying

        public:

        vector<LL>* seq_id_to_color_id;
        ParallelBinaryOutputWriter* out;
        LL output_buffer_max_size;
        LL output_buffer_size;
        char* output_buffer;
        BOSS* boss;
        Coloring* coloring;

        AlignerThread(vector<LL>* seq_id_to_color_id, ParallelBinaryOutputWriter* out, LL output_buffer_max_size, BOSS* boss, Coloring* coloring) 
        : seq_id_to_color_id(seq_id_to_color_id), out(out), output_buffer_max_size(output_buffer_max_size),       output_buffer_size(0), boss(boss), coloring(coloring){
            output_buffer = (char*)malloc(output_buffer_max_size);
        }

        void write(LL node_id, LL color_id){
            LL space_left = output_buffer_max_size - output_buffer_size;

            if(space_left < 8+8){
                // Flush buffer
                out->write(output_buffer, output_buffer_size);
                output_buffer_size = 0;
            }
            
            write_big_endian_LL(output_buffer + output_buffer_size, node_id);
            write_big_endian_LL(output_buffer + output_buffer_size + 8, color_id);
            output_buffer_size += 8+8;
            
        }

        virtual void callback(const char* S, LL S_size, int64_t seq_id){
            LL color = (*seq_id_to_color_id)[seq_id];
            write_log("Adding colors for sequence " + std::to_string(seq_id));
            if(S_size >= boss->get_k()){
                LL node = boss->search(S);
                assert(node >= 0);
                write(node, color);
                for(LL i = boss->get_k(); i < S_size; i++){
                    node = boss->walk(node, S[i]);
                    assert(node >= 0);
                    if(coloring->redundancy_marks[node] == 0){
                        write(node,color);
                    }
                }
            }
        }

        virtual void finish(){
            if(output_buffer_size > 0){
                // Flush buffer
                out->write(output_buffer, output_buffer_size);
                output_buffer_size = 0;
            }
            free(output_buffer);
        }
    };

    // fastafile: the file containing the reference sequences in FASTA format
    // colors_assignments: for each sequence in the fastafile, an integer representing
    // the color of that sequence. The distinct colors should be a contiguous range of
    // integers from 0 to (number of colors -1).
    void add_colors(BOSS& boss, string fastafile, vector<LL> colors_assignments, LL ram_bytes, LL n_threads){

        n_colors = *std::max_element(colors_assignments.begin(), colors_assignments.end()) + 1;
        assert(is_contiguous_integer_range(colors_assignments, 0,n_colors-1));

        write_log("Marking redundant color sets");
        mark_redundant_color_sets(fastafile, boss);

        write_log("Getting (node,color) pairs");
        string node_color_pairs_file = get_node_color_pairs(boss, fastafile, colors_assignments, n_threads);

        write_log("Sorting (node,color) pairs");
        string sorted = EM_sort_big_endian_LL_pairs(node_color_pairs_file, ram_bytes, 0, n_threads);
        temp_file_manager.delete_file(node_color_pairs_file);

        write_log("Deleting duplicate (node,color) pairs");
        string filtered_binary = EM_delete_duplicate_LL_pair_records(sorted);
        temp_file_manager.delete_file(sorted);

        write_log("Collecting color sets");
        string collected = EM_collect_colorsets_binary(filtered_binary);
        temp_file_manager.delete_file(filtered_binary);

        write_log("Sorting color sets");
        string by_colorsets = EM_sort_by_colorsets_binary(collected, ram_bytes, n_threads);
        temp_file_manager.delete_file(collected);

        write_log("Collecting node sets");
        string collected2 = EM_collect_nodes_by_colorset_binary(by_colorsets);
        temp_file_manager.delete_file(by_colorsets);

        write_log("Building packed representation");
        build_packed_representation(collected2, ram_bytes, boss, n_threads);
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

    // Returns a file of pairs (node, color), sorted by node, with possible duplicates
    string get_node_color_pairs(BOSS& boss, string fastafile, vector<LL>& seq_id_to_color_id, LL n_threads){
        //vector<pair<LL,LL> > node_color_pairs;
        string node_color_pair_filename = temp_file_manager.get_temp_file_name("");
        ParallelBinaryOutputWriter writer(node_color_pair_filename);

        vector<DispatcherConsumerCallback*> threads;
        for(LL i = 0; i < n_threads; i++){
            AlignerThread* T = new AlignerThread(&seq_id_to_color_id, &writer, 1024*1024, &boss, this);
            threads.push_back(T);
        }
        run_dispatcher(threads, fastafile, 1024*1024);

        // Clean up
        for(DispatcherConsumerCallback* t : threads) delete t;
        writer.flush();

        return node_color_pair_filename;

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
            if(!redundancy_marks[node])
                class_to_nodes[color_class].push_back(node);
            color_class.clear();
            run_start = run_end + 1;
        }
        return class_to_nodes;
    }

    // Returns pair (number of colorsets, total number of elements in all colorsets)
    pair<LL,LL> count_colorsets(string infile){
        ifstream in(infile, ios::binary);
        LL n_sets = 0;
        LL total_size = 0;
        vector<char> buffer(16);
        while(true){
            in.read(buffer.data(), 16);
            if(!in.good()) break;

            LL record_length = parse_big_endian_LL(buffer.data() + 0);
            LL number_of_nodes = parse_big_endian_LL(buffer.data() + 8);
            LL number_of_colors = (record_length - 16 - number_of_nodes*8)/8;
            n_sets++;
            total_size += number_of_colors;

            // Read the rest (not used)
            while(buffer.size() < record_length)
                buffer.resize(buffer.size() * 2);
            in.read(buffer.data(), record_length - 16); 
        }

        return {n_sets, total_size};
    }

    // Input: pairs (node set, color set)
    void build_packed_representation(string infile, LL ram_bytes, BOSS& boss, LL n_threads){

        // todo: make these private functions return stuff instead of modifying the object state?

        // Concatenate the colorsets and mark the borders
        // Write pairs (node, color set id) to disk.

        LL n_classes, total_size_of_color_sets;
        tie(n_classes, total_size_of_color_sets) = count_colorsets(infile);

        color_sets = sdsl::int_vector<>(total_size_of_color_sets, 0, ceil(log2(n_colors)));
        color_set_starts.resize(color_sets.size());
        sdsl::util::set_to_value(color_set_starts,0);

        nonempty.resize(boss.get_number_of_nodes());
        sdsl::util::set_to_value(nonempty,0);

        nonempty_and_nonredundant.resize(boss.get_number_of_nodes());
        sdsl::util::set_to_value(nonempty_and_nonredundant,0);

        LL index_in_color_sets = 0;
        LL class_id = 0;

        // Iterate all distinct color sets
        //string line;
        ifstream in(infile, ios::binary);
        string node_to_color_id_pairs_filename = temp_file_manager.get_temp_file_name("");
        ofstream node_to_color_id_pairs_out(node_to_color_id_pairs_filename, ios::binary);
        LL n_marks = 0;
        vector<char> buffer(16);

        vector<LL> node_set; // Reusable space
        vector<LL> color_set; // Reusable space
        while(true){

            node_set.clear();
            color_set.clear();

            in.read(buffer.data(), 16);
            if(!in.good()) break;

            LL record_length = parse_big_endian_LL(buffer.data() + 0);
            LL number_of_nodes = parse_big_endian_LL(buffer.data() + 8);
            while(buffer.size() < record_length)
                buffer.resize(buffer.size() * 2);
            in.read(buffer.data() + 16, record_length - 16); // Read the rest

            
            for(LL i = 0; i < number_of_nodes; i++){
                LL node = parse_big_endian_LL(buffer.data() + 16 + i*8);
                node_set.push_back(node);
            }
            
            LL number_of_colors = (record_length - 16 - number_of_nodes*8) / 8;
            for(LL i = 0; i < number_of_colors; i++){
                LL color = parse_big_endian_LL(buffer.data() + 16 + number_of_nodes*8 + i*8);
                color_set.push_back(color);
            }

            // Store the color ids in the color set
            color_set_starts[index_in_color_sets] = 1;
            for(LL x : color_set){
                color_sets[index_in_color_sets] = x;
                index_in_color_sets++;
            }

            // Mark the nodes with this colorset to have a nonempty color set.
            for(LL node : node_set){
                nonempty[node] = 1;
                if(redundancy_marks[node] == 0){ // This should be always true
                    nonempty_and_nonredundant[node] = 1;
                    n_marks++;
                    write_big_endian_LL(node_to_color_id_pairs_out, node);
                    write_big_endian_LL(node_to_color_id_pairs_out, class_id);
                } else{
                    cerr << "Error: this code line should never be reached" << endl; exit(1);
                }
            }
            class_id++;
        }

        node_to_color_id_pairs_out.close();
        
        // Build node_to_color_set_id
        node_to_color_set_id = sdsl::int_vector<>(n_marks, 0, ceil(log2(n_classes)));
        string sorted_out = EM_sort_big_endian_LL_pairs(node_to_color_id_pairs_filename, ram_bytes, 0, n_threads);
        temp_file_manager.delete_file(node_to_color_id_pairs_filename);
        ifstream sorted_in(sorted_out);
        vector<char> buffer2(8+8);
        LL idx = 0;
        while(true){
            sorted_in.read(buffer.data(), 8+8);
            if(!sorted_in.good()) break;
            LL color_set_id = parse_big_endian_LL(buffer.data() + 8);
            node_to_color_set_id[idx] = color_set_id;
            idx++;
        }

        temp_file_manager.delete_file(sorted_out);

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
            tcase.concat += read_separator + tcase.references[i];
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
                        refs.push_back(get_random_dna_string(30, 2));
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
        coloring.add_colors(boss, fastafilename, tcase.seq_id_to_color_id, 1000, 3);

        //cout << tcase.colex_kmers << endl;
        //cout << coloring.to_string_internals() << endl;

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