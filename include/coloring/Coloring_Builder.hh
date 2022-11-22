#pragma once

#include <algorithm>
#include <string>
#include <utility>
#include <vector>
#include <cstdint>
#include <cstring>
#include <variant>
#include <mutex>

#include <sdsl/bit_vectors.hpp>

#include "core_kmer_marker.hh"
#include "backward_traversal.hh"
#include "Sparse_Uint_Array.hh"
#include "Coloring.hh"

#include "sbwt/buffered_streams.hh"
#include "sbwt/EM_sort/bit_level_stuff.hh"
#include "sbwt/EM_sort/EM_sort.hh"
#include "sbwt/globals.hh"
#include "sbwt/SeqIO.hh"

#include "WorkDispatcher.hh"
#include "hybrid_color_set.hh"
#include "Roaring_Color_Set.hh"
#include "Fixed_Width_Int_Color_Set.hh"

#include "ggcat.hh"

using namespace ggcat;

class Colored_Unitig_Stream{ // In-memory implementation for now. Todo: streaming from a ggcat database

    public:

        vector<string> unitigs;
        vector<vector<int64_t> > color_sets;
        int64_t unitig_idx;
        int64_t color_set_idx;

        Colored_Unitig_Stream(const vector<string>& unitigs, const vector<vector<int64_t>>& color_sets):
            unitigs(unitigs), color_sets(color_sets), unitig_idx(0), color_set_idx(0) {
                assert(unitigs.size() == color_sets.size());
        }

        bool done(){
            return unitig_idx >= unitigs.size();
        }

        string next_unitig(){
            return unitigs[unitig_idx++];
        }

        bool next_colors_are_different(){
            return true;
        }

        vector<int64_t> next_colors(){
            return color_sets[color_set_idx++];
        }

};

class Colored_Unitig_Stream_GGCAT{

    vector<string> unitigs;
    vector<vector<int64_t> > color_sets;
    int64_t unitig_idx = 0;
    int64_t color_set_idx = 0;

    public:

        Colored_Unitig_Stream_GGCAT(vector<string>& filenames, int64_t mem_gigas, int64_t n_threads, int64_t k, bool add_rc) {

            GGCATConfig config;

            config.use_temp_dir = true;
            config.temp_dir = get_temp_file_manager().get_dir();
            config.memory = mem_gigas;
            config.prefer_memory = true;
            config.total_threads_count = n_threads;
            config.intermediate_compression_level = -1;

            config.use_stats_file = false;
            config.stats_file = "";

            GGCATInstance *instance = GGCATInstance::create(config);

            std::string graph_file = get_temp_file_manager().create_filename("",".fa");

            std::vector<std::string> color_names;
            for(int64_t i = 0; i < filenames.size(); i++){
                color_names.push_back(to_string(i));
            }

            std::string output_file = instance->build_graph_from_files(
                Slice<std::string>(filenames.data(), filenames.size()),
                graph_file,
                k,
                n_threads,
                !add_rc,
                1,
                ExtraElaborationStep_UnitigLinks,
                true,
                Slice<std::string>(color_names.data(), color_names.size()),
                -1);


            auto file_color_names = GGCATInstance::dump_colors(GGCATInstance::get_colormap_file(graph_file));
            std::mutex print_kmer_lock;

            instance->dump_unitigs(
                graph_file,
                k,
                1,
                // WARNING: this function is called asynchronously from multiple threads, so it must be thread-safe.
                // Also the same_colors boolean is referred to the previous call of this function from the current thread.
                // Number of threads is set to 1 just above, so no lock needed at the moment.
                [&](Slice<char> read, Slice<uint32_t> colors, bool same_colors){
                    std::lock_guard<std::mutex> _lock(print_kmer_lock);
                    try{
                        this->unitigs.push_back(string(read.data, read.data + read.size));
                        this->unitigs.push_back(sbwt::get_rc(unitigs.back())); // Add also the reverse complement

                        if(colors.size == 0){
                            cerr << "BUG: ggcat unitig has empty color set" << endl;
                        }

                        vector<int64_t> colorset;
                        for (size_t i = 0; i < colors.size; i++){
                            colorset.push_back(colors.data[i]);
                        }

                        this->color_sets.push_back(colorset);
                        this->color_sets.push_back(colorset); // Add the same colors for the reverse complement

                    } catch(const std::exception& e){
                        std::cerr << "Caught Error: " << e.what() << '\n';
                        exit(1);
                    }
                    //std::cout << "] same_colors: " << same_colors << std::endl; // TODO
                },
                true);
        }

        bool done(){
            return unitig_idx >= unitigs.size();
        }

        bool next_colors_are_different(){
            return true;
        }

        string next_unitig(){
            return unitigs[unitig_idx++];
        }

        vector<int64_t> next_colors(){
            return color_sets[color_set_idx++];
        }

};

// Color stream from an in-memory vector
class In_Memory_Color_Stream : public Metadata_Stream{

private:

    vector<int64_t> colors;
    int64_t idx = 0; // Current index
    std::array<uint8_t, 8> dummy;

public:

    In_Memory_Color_Stream(vector<int64_t> colors) : colors(colors) {}

    virtual std::array<uint8_t, 8> next(){
        if(idx == colors.size()) return dummy; // Done

        std::array<uint8_t, 8> ret;
        int64_t* ptr = (int64_t*)ret.data(); // Interpret as int64_t
        *ptr = colors[idx++];
        
        return ret;

    }

};


template<typename colorset_t = SDSL_Variant_Color_Set,
         typename sequence_reader_t = sbwt::SeqIO::Reader<>> 
requires Color_Set_Interface<colorset_t>
class Coloring_Builder{

private:

    class ColorPairAlignerThread : public DispatcherConsumerCallback {
        ParallelBinaryOutputWriter& out;
        const std::size_t output_buffer_max_size;
        std::size_t output_buffer_size = 0;
        char* output_buffer;
        const plain_matrix_sbwt_t& index;
        const sdsl::bit_vector& cores;
        std::int64_t largest_color_id = 0;

    public:
        ColorPairAlignerThread(ParallelBinaryOutputWriter& out,
                      const std::size_t output_buffer_max_size,
                      const plain_matrix_sbwt_t& index,
                      const sdsl::bit_vector& cores) :
            out(out),
            output_buffer_max_size(output_buffer_max_size),
            index(index),
            cores(cores) {
            output_buffer = new char[output_buffer_max_size];
        }

        ColorPairAlignerThread(const ColorPairAlignerThread&) = delete;
        ColorPairAlignerThread& operator=(const ColorPairAlignerThread&) = delete;

        void write(const std::int64_t node_id, const std::int64_t color_id) {
            const std::size_t space_left = output_buffer_max_size - output_buffer_size;

            if (space_left < 8+8) {
                out.write(output_buffer, output_buffer_size);
                output_buffer_size = 0;
            }

            write_big_endian_LL(output_buffer + output_buffer_size, node_id);
            write_big_endian_LL(output_buffer + output_buffer_size + 8, color_id);
            output_buffer_size += 8+8;
        }

        // It is our responsibility to interpret AND DELETE the metadata
        virtual void callback(const char* S,
                              int64_t S_size,
                              int64_t string_id,
                              std::array<uint8_t, 8> metadata) {

            int64_t color = *reinterpret_cast<int64_t*>(metadata.data()); // Interpret as int64_t
            const std::size_t k = index.get_k();

            write_log("Adding colors for sequence " + std::to_string(string_id), LogLevel::MINOR);
            if (S_size >= k) {
                const auto res = index.streaming_search(S, S_size);
                for (const auto node : res) {
                    if (node >= 0 && cores[node] == 1) {
                        write(node, color);
                        largest_color_id = std::max(largest_color_id, color);
                    }
                }
            }
        }

        virtual void finish() {
            if (output_buffer_size > 0) {
                out.write(output_buffer, output_buffer_size);
                output_buffer_size = 0;
            }

            delete[] output_buffer;
        }

        std::int64_t get_largest_color_id() const {
            return largest_color_id;
        }
    };

    // Return the filename of the generated node-color pairs, and the largest color id
    pair<std::string, int64_t> get_node_color_pairs(const plain_matrix_sbwt_t& index,
                                     sequence_reader_t& reader,
                                     Metadata_Stream* metadata_stream,
                                     const sdsl::bit_vector& cores,
                                     const std::size_t n_threads) {
        const std::string outfile = get_temp_file_manager().create_filename();

        ParallelBinaryOutputWriter writer(outfile);

        std::vector<DispatcherConsumerCallback*> threads;
        for (std::size_t i = 0; i < n_threads; ++i) {
            ColorPairAlignerThread* T = new ColorPairAlignerThread(
                                                 writer,
                                                 1024*1024,
                                                 index,
                                                 cores);
            threads.push_back(T);
        }

        run_dispatcher(threads, reader, metadata_stream, 1024*1024);

        std::vector<std::int64_t> largest_color_ids;
        for (DispatcherConsumerCallback* t : threads) {
            ColorPairAlignerThread* cpat = static_cast<ColorPairAlignerThread*>(t);
            largest_color_ids.push_back(cpat->get_largest_color_id());
            delete t;
        }

        int64_t largest_color_id = *std::max_element(largest_color_ids.begin(), largest_color_ids.end());

        writer.flush();

        return {outfile, largest_color_id};
    }

    std::string delete_duplicate_pairs(const std::string& infile) {
        std::string outfile = get_temp_file_manager().create_filename();

        Buffered_ifstream<> in(infile, std::ios::binary);
        Buffered_ofstream<> out(outfile, std::ios::binary);

        char prev[8+8]; // two long longs
        char cur[8+8]; // two long longs

        std::size_t record_count = 0;
        while (in.read(cur, 8+8)){

            if (record_count == 0 || std::memcmp(prev, cur, 8+8) != 0){
                // The first record or different from the previous record
                out.write(cur, 8+8);
            }

            std::memcpy(prev, cur, 8+8);
            ++record_count;
        }
        out.flush();

        return outfile;
    }

    std::string collect_colorsets(const std::string& infile){
        std::string outfile = get_temp_file_manager().create_filename();

        Buffered_ifstream<> in(infile, ios::binary);
        Buffered_ofstream<> out(outfile, ios::binary);

        std::int64_t active_key = -1;
        std::vector<std::int64_t> cur_value_list;

        char buffer[8+8];

        while (true) {
            in.read(buffer,8+8);
            if (in.eof()) break;

            std::int64_t key = parse_big_endian_LL(buffer + 0);
            std::int64_t value = parse_big_endian_LL(buffer + 8);

            if (key == active_key)
                cur_value_list.push_back(value);
            else {
                if (active_key != -1) {
                    std::sort(cur_value_list.begin(), cur_value_list.end());
                    std::int64_t record_size = 8 * (1 + 1 + cur_value_list.size());
                    write_big_endian_LL(out, record_size);
                    write_big_endian_LL(out, active_key);

                    for (auto x : cur_value_list){
                        write_big_endian_LL(out, x);
                    }
                }

                active_key = key;
                cur_value_list.clear();
                cur_value_list.push_back(value);
            }
        }

        // Last one
        if(active_key != -1){
            sort(cur_value_list.begin(), cur_value_list.end());
            std::int64_t record_size = 8 * (1 + 1 + cur_value_list.size());
            write_big_endian_LL(out, record_size);
            write_big_endian_LL(out, active_key);

            for (auto x : cur_value_list){
                write_big_endian_LL(out, x);
            }
        }

        out.flush();

        return outfile;
    }

    std::string sort_by_colorsets(const std::string& infile,
                                  const std::int64_t ram_bytes,
                                  const std::int64_t n_threads) {
        auto cmp = [&](const char* x, const char* y) -> bool{
            std::int64_t nx = parse_big_endian_LL(x);
            std::int64_t ny = parse_big_endian_LL(y);
            std::int64_t c = std::memcmp(x + 8 + 8, y + 8 + 8, std::min(nx-8-8,ny-8-8));

            if (c < 0)
                return true;
            else if (c > 0)
                return false;
            else { // c == 0
                return nx < ny;
            }
        };

        std::string outfile = get_temp_file_manager().create_filename();
        EM_sort_variable_length_records(infile, outfile, cmp, ram_bytes, n_threads);
        return outfile;
    }

    std::string collect_nodes_by_colorset(const std::string& infile){
        std::string outfile = get_temp_file_manager().create_filename();

        Buffered_ifstream<> in(infile, ios::binary);
        Buffered_ofstream<> out(outfile, ios::binary);

        std::vector<std::int64_t> active_key;
        std::vector<std::int64_t> cur_value_list;

        std::vector<char> buffer(8);
        std::vector<std::int64_t> key; // Reusable space
        while (true) {
            key.clear();
            in.read(buffer.data(),8);

            if (in.eof())
                break;

            std::int64_t record_len = parse_big_endian_LL(buffer.data());

            assert(record_len >= 8+8+8);
            while (buffer.size() < record_len)
                buffer.resize(buffer.size()*2);

            in.read(buffer.data()+8,record_len-8); // Read the rest
            std::int64_t value = parse_big_endian_LL(buffer.data()+8);
            std::int64_t color_set_size = (record_len - 8 - 8) / 8;

            for (std::int64_t i = 0; i < color_set_size; i++) {
                key.push_back(parse_big_endian_LL(buffer.data() + 8 + 8 + i*8));
            }

            if (key == active_key)
                cur_value_list.push_back(value);
            else {
                if (active_key.size() != 0) {
                    std::sort(cur_value_list.begin(), cur_value_list.end());
                    std::int64_t record_size = 8 + 8 + 8*cur_value_list.size() + 8*active_key.size();
                    assert(record_size > 0);

                    // record = (record length, number of nodes, node list, color list)
                    write_big_endian_LL(out, record_size);
                    write_big_endian_LL(out, cur_value_list.size());

                    for (auto x : cur_value_list) {
                        write_big_endian_LL(out, x);
                    }

                    for (auto x : active_key) {
                        write_big_endian_LL(out, x);
                    }
                }

                active_key = key;
                cur_value_list.clear();
                cur_value_list.push_back(value);
            }
        }

        // Last one
        if (active_key.size() != 0) {
            std::sort(cur_value_list.begin(), cur_value_list.end());
            std::int64_t record_size = 8 + 8 + 8*cur_value_list.size() + 8*active_key.size();

            write_big_endian_LL(out, record_size);
            write_big_endian_LL(out, cur_value_list.size());

            for (auto x : cur_value_list) {
                write_big_endian_LL(out, x);
            }

            for (auto x : active_key) {
                write_big_endian_LL(out, x);
            }
        }

        out.flush();

        return outfile;
    }

    void build_representation(Coloring<colorset_t>& coloring, const std::string& infile, const sdsl::bit_vector& cores, int64_t colorset_sampling_distance, int64_t ram_bytes, int64_t n_threads) {

        SBWT_backward_traversal_support backward_support(coloring.index_ptr);

        Buffered_ifstream<> in(infile, ios::binary);
        vector<char> buffer(16);

        vector<std::int64_t> node_set; // Reusable space
        vector<std::int64_t> colors_set; // Reusable space

        std::size_t set_id = 0;
        Sparse_Uint_Array_Builder builder(cores.size(), ram_bytes, n_threads);

        auto callback = [&](int64_t node){
            builder.add(node, set_id);;
        };

        while (true) {
            node_set.clear();
            colors_set.clear();

            in.read(buffer.data(), 16);
            if (in.eof())
                break;

            const auto record_length = parse_big_endian_LL(buffer.data() + 0);
            const auto number_of_nodes = parse_big_endian_LL(buffer.data() + 8);
            while (buffer.size() < record_length)
                buffer.resize(buffer.size() * 2);
            in.read(buffer.data() + 16, record_length - 16); // Read the rest

            for (std::int64_t i = 0; i < number_of_nodes; i++) {
                std::int64_t node = parse_big_endian_LL(buffer.data() + 16 + i*8);
                node_set.push_back(node);
            }

            std::int64_t number_of_colors = (record_length - 16 - number_of_nodes*8) / 8;
            for (std::int64_t i = 0; i < number_of_colors; i++) {
                std::int64_t color = parse_big_endian_LL(buffer.data() + 16 + number_of_nodes*8 + i*8);
                colors_set.push_back(color);
            }

            coloring.sets.add_set(colors_set);
            coloring.total_color_set_length += colors_set.size();

            for (int64_t node : node_set) {
                builder.add(node, set_id);
                iterate_unitig_node_samples(cores, backward_support, node, colorset_sampling_distance, callback);
            }

            ++set_id;
        }

        coloring.node_id_to_color_set_id = builder.finish();
        coloring.sets.prepare_for_queries();

    }

    // Walks backward from from_node and marks every colorset_sampling_distance node on the way
    // Calls the given callback at every sampled node. The callback takes the node id of the sampled node.
    template<typename callback_t>
    void iterate_unitig_node_samples(const sdsl::bit_vector& cores, SBWT_backward_traversal_support& backward_support, int64_t from_node, int64_t colorset_sampling_distance, callback_t&& callback){
        assert(cores[from_node] == 1);
        int64_t in_neighbors[4];
        int64_t indegree;
        backward_support.list_DBG_in_neighbors(from_node, in_neighbors, indegree);
        for(int64_t i = 0; i < indegree; i++){
            int64_t u = in_neighbors[i];
            int64_t counter = 0;
            while(cores[u] == 0){
                counter++;
                if(counter == colorset_sampling_distance){
                    callback(u);
                    counter = 0;
                }
                backward_support.list_DBG_in_neighbors(u, in_neighbors, indegree);
                if(indegree == 0) break; // Root node
                if(indegree >= 2) break; // Predecessors are already marked
                u = in_neighbors[0]; // The only in-neighbor
                // This loop is guaranteed to terminate. Proof:
                // If the in-degree and out-degree of u always stays 1, then we come back to the original node eventually,
                // which is marked, so we stop. Otherwise, the in-degree becomes 0 or >= 2 at some point, and we hit
                // the break statements above.
            }
        }
    }

    public:


    void build_coloring(
                    Coloring<colorset_t>& coloring,
                    const plain_matrix_sbwt_t& index,
                    sequence_reader_t& sequence_reader,
                    const vector<int64_t>& color_assignment,
                    const std::int64_t ram_bytes,
                    const std::int64_t n_threads,
                    int64_t colorset_sampling_distance) {
        In_Memory_Color_Stream imcs(color_assignment);
        build_coloring(coloring, index, sequence_reader, &imcs, ram_bytes, n_threads, colorset_sampling_distance);
    }

    void build_coloring(
                    Coloring<colorset_t>& coloring,
                    const plain_matrix_sbwt_t& index,
                    sequence_reader_t& sequence_reader,
                    Metadata_Stream* metadata_stream,
                    const std::int64_t ram_bytes,
                    const std::int64_t n_threads,
                    int64_t colorset_sampling_distance) {

        coloring.index_ptr = &index;

        write_log("Marking core kmers", LogLevel::MAJOR);
        core_kmer_marker<sequence_reader_t> ckm;
        ckm.mark_core_kmers(sequence_reader, index);
        sdsl::bit_vector cores = ckm.core_kmer_marks;

        sequence_reader.rewind_to_start(); // Need this reader again for node-colors pairs

        write_log("Getting node color pairs", LogLevel::MAJOR);
        std::string node_color_pairs; int64_t largest_color_id;
        std::tie(node_color_pairs, largest_color_id) = get_node_color_pairs(index, sequence_reader, metadata_stream, cores, n_threads);
        coloring.largest_color_id = largest_color_id;

        write_log("Sorting node color pairs", LogLevel::MAJOR);
        const std::string sorted_pairs = get_temp_file_manager().create_filename();

        auto cmp = [&](const char* A, const char* B) -> bool {
            std::int64_t x_1, y_1, x_2, y_2;
            x_1 = parse_big_endian_LL(A + 0);
            y_1 = parse_big_endian_LL(A + 8);
            x_2 = parse_big_endian_LL(B + 0);
            y_2 = parse_big_endian_LL(B + 8);

            return std::make_pair(x_1, y_1) < std::make_pair(x_2, y_2);
        };

        EM_sort_constant_binary(node_color_pairs, sorted_pairs, cmp, ram_bytes, 16, n_threads);
        get_temp_file_manager().delete_file(node_color_pairs);

        write_log("Removing duplicate node color pairs", LogLevel::MAJOR);
        const std::string filtered_pairs = delete_duplicate_pairs(sorted_pairs);
        get_temp_file_manager().delete_file(sorted_pairs);

        write_log("Collecting colors", LogLevel::MAJOR);
        const std::string collected_sets = collect_colorsets(filtered_pairs);
        get_temp_file_manager().delete_file(filtered_pairs);

        write_log("Sorting color sets", LogLevel::MAJOR);
        const std::string sorted_sets = sort_by_colorsets(collected_sets, ram_bytes, n_threads);
        get_temp_file_manager().delete_file(collected_sets);

        write_log("Collecting nodes", LogLevel::MAJOR);
        string collected_nodes = collect_nodes_by_colorset(sorted_sets);
        get_temp_file_manager().delete_file(sorted_sets);

        write_log("Building representation", LogLevel::MAJOR);
        build_representation(coloring, collected_nodes, cores, colorset_sampling_distance, ram_bytes, n_threads);
        get_temp_file_manager().delete_file(collected_nodes);

        write_log("Representation built", LogLevel::MAJOR);
    }

    template<typename colored_unitig_stream_t>
    void build_from_colored_unitigs(Coloring<colorset_t>& coloring,
                    sequence_reader_t& sequence_reader, // The original sequences, not the unitigs. Used to mark core k-mers
                    const plain_matrix_sbwt_t& SBWT,
                    const std::int64_t ram_bytes,
                    const std::int64_t n_threads,
                    int64_t colorset_sampling_distance,
                    colored_unitig_stream_t& colored_unitig_stream){
        
        coloring.index_ptr = &SBWT;

        write_log("Marking core kmers", LogLevel::MAJOR);
        core_kmer_marker<sequence_reader_t> ckm;
        ckm.mark_core_kmers(sequence_reader, SBWT);
        sdsl::bit_vector cores = ckm.core_kmer_marks;

        SBWT_backward_traversal_support backward_support(coloring.index_ptr);

        Sparse_Uint_Array_Builder builder(cores.size(), ram_bytes, n_threads);

        int64_t color_set_id = -1;

        auto add_color_set_pointer = [&](int64_t node_id){
            builder.add(node_id, color_set_id);
        };

        coloring.largest_color_id = -1;
        vector<int64_t> colors;
        while(!colored_unitig_stream.done()){
            string unitig = colored_unitig_stream.next_unitig();
            if(colored_unitig_stream.next_colors_are_different()){
                // Read new color set
                colors = colored_unitig_stream.next_colors();
                for(int64_t x : colors)
                    coloring.largest_color_id = max(x, coloring.largest_color_id);
                
                color_set_id++;
            }

            
            // Store color set
            coloring.sets.add_set(colors);
            coloring.total_color_set_length += colors.size();

            // Store pointers to the color set
            for(int64_t colex_rank : SBWT.streaming_search(unitig)){
                if(colex_rank != -1 && cores[colex_rank]){
                    add_color_set_pointer(colex_rank);
                    iterate_unitig_node_samples(cores, backward_support, colex_rank, colorset_sampling_distance, add_color_set_pointer);
                }
            }
        }

        coloring.node_id_to_color_set_id = builder.finish();
        coloring.sets.prepare_for_queries();        
    }
};
