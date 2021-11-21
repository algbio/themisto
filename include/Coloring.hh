#pragma once

#include "libwheeler/BOSS.hh"
#include "TempFileManager.hh"
#include "EM_algorithms.hh"
#include "WorkDispatcher.hh"
#include "libwheeler/BOSS_builder.hh"
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
// Every color name is assigned a distinct id, starting from
// zero such that the new ids are handed out in the order of 
// appearance of the names in the user-provided colorfile.

// The input for indexing is:
// - A fasta file with m sequences
// - A colorfile with m lines such that line i gives the
//   name of the color of the i-th sequence

// The input for querying is:
// - A fasta file of reads 
// The output for querying is:
// - For each read, a set of color names (not color ids)

private:

    // private copy construtor and assignment operator: no copying, 
    // else might invalidate the pointer in the rank/select supports
    Coloring(const Coloring&); 
    Coloring& operator=(const Coloring&); 

public:

    inline static const LL VERSION = 1; // Increment this after non-breaking changes
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
        BOSS<sdsl::bit_vector>* boss;
        Coloring* coloring;

        AlignerThread(vector<LL>* seq_id_to_color_id, ParallelBinaryOutputWriter* out, LL output_buffer_max_size, BOSS<sdsl::bit_vector>* boss, Coloring* coloring) 
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
            if(seq_id >= seq_id_to_color_id->size()){
                throw std::runtime_error("Error at: seq_id >= seq_id_to_color_id->size()");
            }
            LL color = (*seq_id_to_color_id)[seq_id];
            write_log("Adding colors for sequence " + std::to_string(seq_id));
            if(S_size >= boss->get_k() + 1){
                // We have a +1 in the if condition because the graph is edge centric (= defined by (k+1)-mers),
                // so sequences of size exactly k should not produce colors.
                LL node = boss->find_kmer(S);
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
    void add_colors(BOSS<sdsl::bit_vector>& boss, string fastafile, vector<LL> colors_assignments, LL ram_bytes, LL n_threads, LL colorset_sampling_distance){

        n_colors = *std::max_element(colors_assignments.begin(), colors_assignments.end()) + 1;

        write_log("Marking redundant color sets");
        mark_redundant_color_sets(fastafile, boss);

        write_log("Getting (node,color) pairs");
        string node_color_pairs_file = get_node_color_pairs(boss, fastafile, colors_assignments, n_threads);

        write_log("Sorting (node,color) pairs");
        string sorted = EM_sort_big_endian_LL_pairs(node_color_pairs_file, ram_bytes, 0, n_threads);
        get_temp_file_manager().delete_file(node_color_pairs_file);

        write_log("Deleting duplicate (node,color) pairs");
        string filtered_binary = EM_delete_duplicate_LL_pair_records(sorted);
        get_temp_file_manager().delete_file(sorted);

        write_log("Collecting color sets");
        string collected = EM_collect_colorsets_binary(filtered_binary);
        get_temp_file_manager().delete_file(filtered_binary);

        write_log("Sorting color sets");
        string by_colorsets = EM_sort_by_colorsets_binary(collected, ram_bytes, n_threads);
        get_temp_file_manager().delete_file(collected);

        write_log("Collecting node sets");
        string collected2 = EM_collect_nodes_by_colorset_binary(by_colorsets);
        get_temp_file_manager().delete_file(by_colorsets);

        write_log("Building packed representation");
        build_packed_representation(collected2, ram_bytes, boss, colorset_sampling_distance, n_threads);
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

    string to_string(BOSS<sdsl::bit_vector>& boss){
        stringstream ss;
        for(LL node = 0; node < nonempty.size(); node++){
            ss << node << ": ";
            for(LL color : get_colorset(node, boss)) ss << color << " ";
            ss << "\n";
        }
        
        return ss.str();
    }

    // Buffer must have enough space to accommodate all colors.
    // Puts the colors into the buffer starting from index 0 and
    // returns the number of colors added.
    // todo: tests
    LL get_colorset_to_buffer(LL node, BOSS<sdsl::bit_vector>& boss, vector<LL>& buffer){
        LL id = get_colorset_id(node, boss);
        return get_colorset_to_buffer_by_id(id, buffer);
    }

    // Returns the number of elements put into the buffer. Buffer should be big enough to accommodate
    // all possible distinct colors.
    LL get_colorset_to_buffer_by_id(LL colorset_id, vector<LL>& buffer){
        if(colorset_id == -1) return 0;

        LL start = color_set_starts_ss.select(colorset_id+1);
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

    // Returns -1 if the colorset is empty. Otherwise returns the id of the colorset.
    LL get_colorset_id(LL node, BOSS<sdsl::bit_vector>& boss){
        if(node == 0) return -1; // Root
        while(redundancy_marks[node] == 1){
            node = boss.edge_source(boss.inedge_range(node).first); // Go to next non-redundant
            if(node == 0) return -1; // Root
        }

        // The nodes that have an explicit color set are those that have a non-empty
        // colorset and are non-redundant.
        if(nonempty[node] == 0) {
            return -1;
        } else{
            LL rank = nonempty_and_nonredundant_rs.rank(node);
            return node_to_color_set_id[rank];
        }
    }

    // Super inefficient because allocates a new buffer every time
    set<LL> get_colorset(LL node, BOSS<sdsl::bit_vector>& boss){
        vector<LL> buffer(n_colors);
        LL count = get_colorset_to_buffer(node, boss, buffer);
        set<LL> result(buffer.begin(), buffer.begin() + count);
        return result;
    }

    // Super inefficient because allocates a new buffer every time
    vector<LL> get_colorvec(LL node, BOSS<sdsl::bit_vector>& boss){
        vector<LL> buffer(n_colors);
        LL count = get_colorset_to_buffer(node, boss, buffer);
        buffer.resize(count);
        return buffer;
    }

    void load(istream& is){
        LL stored_version_number;
        is.read((char*)&stored_version_number, sizeof(stored_version_number));
        if(stored_version_number != VERSION)
            throw std::runtime_error("The coloring data structure was built with an incompatible version of the software");

        color_sets.load(is);
        node_to_color_set_id.load(is);
        nonempty.load(is);
        nonempty_rs.load(is);
        color_set_starts.load(is);
        color_set_starts_ss.load(is);
        redundancy_marks.load(is);
        nonempty_and_nonredundant.load(is);
        nonempty_and_nonredundant_rs.load(is);

        color_set_starts_ss.set_vector(&color_set_starts);
        nonempty_rs.set_vector(&nonempty);
        nonempty_and_nonredundant_rs.set_vector(&nonempty_and_nonredundant);

        n_colors = 0;
        for(LL i = 0; i < color_sets.size(); i++) n_colors = max(n_colors, (LL)(color_sets[i]+1));
    }

    LL serialize(ostream& os){
        LL written = 0;

        os.write((char*)&VERSION, sizeof(VERSION));
        written += sizeof(VERSION);

        written += color_sets.serialize(os);
        written += node_to_color_set_id.serialize(os);
        written += nonempty.serialize(os);
        written += nonempty_rs.serialize(os);
        written += color_set_starts.serialize(os);
        written += color_set_starts_ss.serialize(os);
        written += redundancy_marks.serialize(os);
        written += nonempty_and_nonredundant.serialize(os);
        written += nonempty_and_nonredundant_rs.serialize(os);
        return written;
    }

    LL get_number_of_distinct_colorsets(){
        LL ans = 0;
        for(LL i = 0; i < color_set_starts.size(); i++)
            ans += color_set_starts[i];
        return ans;
    }

    LL get_color_set_concatenation_length(){
        return color_set_starts.size();
    }


private:

    bool is_contiguous_integer_range(vector<LL>& v, LL from, LL to){
        set<LL> S(v.begin(), v.end());
        return *S.begin() >= from && *S.rbegin() <= to && S.size() == to - from + 1;
    }

    // Returns a file of pairs (node, color), sorted by node, with possible duplicates
    string get_node_color_pairs(BOSS<sdsl::bit_vector>& boss, string fastafile, vector<LL>& seq_id_to_color_id, LL n_threads){
        //vector<pair<LL,LL> > node_color_pairs;
        string node_color_pair_filename = get_temp_file_manager().create_filename();
        ParallelBinaryOutputWriter writer(node_color_pair_filename);

        vector<DispatcherConsumerCallback*> threads;
        for(LL i = 0; i < n_threads; i++){
            AlignerThread* T = new AlignerThread(&seq_id_to_color_id, &writer, 1024*1024, &boss, this);
            threads.push_back(T);
        }

        Sequence_Reader_Buffered sr(fastafile, FASTA_MODE);
        run_dispatcher(threads, sr, 1024*1024);

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
        Buffered_ifstream in(infile, ios::binary);
        LL n_sets = 0;
        LL total_size = 0;
        vector<char> buffer(16);
        while(true){
            in.read(buffer.data(), 16);
            if(in.eof()) break;

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
    void build_packed_representation(string infile, LL ram_bytes, BOSS<sdsl::bit_vector>& boss, LL colorset_sampling_distance, LL n_threads){

        // todo: make these private functions return stuff instead of modifying the object state?

        // Concatenate the colorsets and mark the borders
        // Write pairs (node, color set id) to disk.

        LL n_classes, total_size_of_color_sets;
        tie(n_classes, total_size_of_color_sets) = count_colorsets(infile);

        color_sets = sdsl::int_vector<>(total_size_of_color_sets, 0, ceil(log2(n_colors)));
        color_set_starts.resize(color_sets.size());
        sdsl::util::set_to_value(color_set_starts,0);

        nonempty.resize(boss.number_of_nodes());
        sdsl::util::set_to_value(nonempty,0);

        nonempty_and_nonredundant.resize(boss.number_of_nodes());
        sdsl::util::set_to_value(nonempty_and_nonredundant,0);

        // We have a temporary vector for marking new nodes. Otherwise we would mark in the same vector
        // that we are using to detect marking path starting points, which would mess things up.
        sdsl::bit_vector redundancy_marks_temp = redundancy_marks;

        LL index_in_color_sets = 0;
        LL class_id = 0;

        // Iterate all distinct color sets
        //string line;
        Buffered_ifstream in(infile, ios::binary);
        string node_to_color_id_pairs_filename = get_temp_file_manager().create_filename();
        Buffered_ofstream node_to_color_id_pairs_out(node_to_color_id_pairs_filename, ios::binary);
        LL n_marks = 0;
        vector<char> buffer(16);

        vector<LL> node_set; // Reusable space
        vector<LL> color_set; // Reusable space
        while(true){

            node_set.clear();
            color_set.clear();

            in.read(buffer.data(), 16);
            if(in.eof()) break;

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
                assert(redundancy_marks[node] == 0);
                if(redundancy_marks[node] == 0){
                    nonempty_and_nonredundant[node] = 1;
                    n_marks++;
                    write_big_endian_LL(node_to_color_id_pairs_out, node);
                    write_big_endian_LL(node_to_color_id_pairs_out, class_id);
                    n_marks += propagate_colorset_pointers(redundancy_marks_temp, node, class_id, colorset_sampling_distance, boss, node_to_color_id_pairs_out);
                } else{
                    cerr << "Error: this code line should never be reached" << endl; exit(1);
                }
            }
            class_id++;
        }

        node_to_color_id_pairs_out.close();
        
        redundancy_marks = redundancy_marks_temp;
        
        // Build node_to_color_set_id
        node_to_color_set_id = sdsl::int_vector<>(n_marks, 0, ceil(log2(n_classes)));
        string sorted_out = EM_sort_big_endian_LL_pairs(node_to_color_id_pairs_filename, ram_bytes, 0, n_threads);
        get_temp_file_manager().delete_file(node_to_color_id_pairs_filename);
        Buffered_ifstream sorted_in(sorted_out);
        vector<char> buffer2(8+8);
        LL idx = 0;
        while(true){
            sorted_in.read(buffer.data(), 8+8);
            if(sorted_in.eof()) break;
            LL color_set_id = parse_big_endian_LL(buffer.data() + 8);
            node_to_color_set_id[idx] = color_set_id;
            idx++;
        }

        get_temp_file_manager().delete_file(sorted_out);

        sdsl::util::init_support(nonempty_rs, &nonempty);
        sdsl::util::init_support(color_set_starts_ss, &color_set_starts);
        sdsl::util::init_support(nonempty_and_nonredundant_rs, &nonempty_and_nonredundant);
    }

    void mark_redundant_color_sets(string fastafile, BOSS<sdsl::bit_vector>& boss){
        // Mark all nodes that fulfill at least one of the following
        // (1) the node is first k-mer of an input sequence
        // (2) a predecessor of the node is the last k-mer of a reference sequence
        // (3) have two or more incoming edges
        // (4) have a predecessor that has two or more outgoing edges
        redundancy_marks.resize(boss.number_of_nodes());
        sdsl::util::set_to_value(redundancy_marks,1); // 1 = redundant, 0 = non-redundant.

        LL k = boss.get_k();

        // Handle cases (1) and (2)
        Sequence_Reader_Buffered sr(fastafile, FASTA_MODE);
        while(true){
            LL len = sr.get_next_read_to_buffer();
            if(len == 0) break;
            if(len >= k+1){ // Reads shorter than k+1 do not contribute to the graph
                // condition (1)
                LL first_node = boss.find_kmer(sr.read_buf); // Guaranteed to be found in the graph
                redundancy_marks[first_node] = 0;

                // condition (2)
                LL last_node = boss.find_kmer(sr.read_buf + len - k);
                pair<int64_t, int64_t> I = boss.outlabel_range(last_node);
                for(LL i = I.first; i <= I.second; i++){
                    LL next = boss.edge_destination(boss.outedge_index_to_wheeler_rank(i));
                    redundancy_marks[next] = 0;
                }
            }
        }

        // Handle cases (3) and (4)
        for(LL v = 1; v < boss.number_of_nodes(); v++){
            // Starting from 1: don't consider the root node 0. This way we always have a predecessor
            pair<int64_t, int64_t> inedge_range = boss.inedge_range(v);
            if(inedge_range.second > inedge_range.first){ // Indegree >= 2
                redundancy_marks[v] = 0;
                continue;
            } else{
                assert(inedge_range.first == inedge_range.second);
                LL u = boss.edge_source(inedge_range.first);
                if(boss.outdegree(u) > 1) redundancy_marks[v] = 0;
            }
        }
    }


    // Returns the number of new marks
    LL propagate_colorset_pointers(sdsl::bit_vector& redundancy_marks_temp, LL from_node, LL class_id, LL colorset_sampling_distance, BOSS<sdsl::bit_vector>& boss, Buffered_ofstream& out){

        LL n_new_marks = 0;
        
        assert(redundancy_marks[from_node] == 0);
        pair<int64_t, int64_t> I = boss.outlabel_range(from_node);
        for(LL i = I.first; i <= I.second; i++){
            LL u = boss.edge_destination(boss.outedge_index_to_wheeler_rank(i));
            LL counter = 0;
            while(redundancy_marks[u] == 1){
                counter++;
                if(counter == colorset_sampling_distance){
                    redundancy_marks_temp[u] = 0; // new mark
                    nonempty_and_nonredundant[u] = 1;
                    nonempty[u] = 1;
                    counter = 0;
                    n_new_marks++;
                    write_big_endian_LL(out, u);
                    write_big_endian_LL(out, class_id);
                }
                pair<int64_t, int64_t> I_u = boss.outlabel_range(u);
                if(I_u == make_pair<int64_t, int64_t>(1,0)) break; // Outdegree is zero
                u = boss.edge_destination(boss.outedge_index_to_wheeler_rank(I_u.first));
                // This loop is guaranteed to terminate. Proof:
                // If the in-degree and out-degree of u always stays 1, then we come back to v eventually,
                // which is marked as non-redundant, so we stop. Otherwise:
                // - If the out-degree of u becomes greater than 1 at some point, then the
                // successors of that u will be marked non-redundant by case (4)
                // - If the in-degree of u becomes greater than 1 at some point, then it is
                //   marked as non-redundant by case (3)
            }
        }
        
        return n_new_marks;
    }

};
