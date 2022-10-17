#pragma once

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include <cstdint>
#include <cstring>

#include <sdsl/bit_vectors.hpp>

#include "core_kmer_marker.hh"
#include "backward_traversal.hh"

#include "Sparse_Uint_Array.hh"
#include "sbwt/buffered_streams.hh"
#include "sbwt/EM_sort/bit_level_stuff.hh"
#include "sbwt/EM_sort/EM_sort.hh"
#include "sbwt/globals.hh"
#include "sbwt/SeqIO.hh"

#include <vector>

#include <cstdint>

#include "WorkDispatcher.hh"
#include "sdsl_color_set.hh"
#include "Roaring_Color_Set.hh"
#include <variant>

template <typename T>
concept Color_Set_Interface = requires(T& t, std::ostream& os, std::istream& is)
{
    { t.size() } -> std::same_as<int64_t>;
    { t.size_in_bits() } -> std::same_as<int64_t>;
    { t.contains(int64_t()) } -> std::same_as<bool>;
    { t.intersection(t) } -> std::same_as<T>;
    { t.do_union(t) } -> std::same_as<T>;
    { t.serialize(os) } -> std::same_as<int64_t>;
    { t.load(is) } -> std::same_as<void>;
    { t.get_colors_as_vector() } -> std::same_as<std::vector<int64_t>>;
};


// Takes as parameter a class that encodes a single color set
template<typename colorset_t = Bitmap_Or_Deltas_ColorSet> requires Color_Set_Interface<colorset_t>
class Coloring {

public:
    class WrongTemplateParameterException : public std::exception{
        const char * what() const noexcept override{
            return "Template type id in a serialized Coloring structure does not match the class template parameter.";
        }
    };

    typedef colorset_t colorset_type;

private:

    std::vector<colorset_t> sets;

    Sparse_Uint_Array node_id_to_color_set_id;
    const plain_matrix_sbwt_t* index_ptr;
    int64_t largest_color_id = 0;
    int64_t total_color_set_length = 0;

    class ColorPairAlignerThread : public DispatcherConsumerCallback {
        const std::vector<std::int64_t>& seq_id_to_color_id;
        ParallelBinaryOutputWriter& out;
        const std::size_t output_buffer_max_size;
        std::size_t output_buffer_size = 0;
        char* output_buffer;
        const plain_matrix_sbwt_t& index;
        const Coloring& coloring;
        const sdsl::bit_vector& cores;
        std::int64_t largest_color_id = 0;

    public:
        ColorPairAlignerThread(const std::vector<std::int64_t>& seq_id_to_color_id,
                      ParallelBinaryOutputWriter& out,
                      const std::size_t output_buffer_max_size,
                      const plain_matrix_sbwt_t& index,
                      const Coloring& coloring,
                      const sdsl::bit_vector& cores) :
            seq_id_to_color_id(seq_id_to_color_id),
            out(out),
            output_buffer_max_size(output_buffer_max_size),
            index(index),
            coloring(coloring),
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

        virtual void callback(const char* S,
                              LL S_size,
                              int64_t string_id) {
            const std::int64_t color = seq_id_to_color_id.at(string_id);
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

    std::string get_node_color_pairs(const plain_matrix_sbwt_t& index,
                                     const std::string& fasta_file,
                                     const std::vector<std::int64_t>& seq_id_to_color_id,
                                     const sdsl::bit_vector& cores,
                                     const std::size_t n_threads) {
        const std::string outfile = get_temp_file_manager().create_filename();

        ParallelBinaryOutputWriter writer(outfile);

        std::vector<DispatcherConsumerCallback*> threads;
        for (std::size_t i = 0; i < n_threads; ++i) {
            ColorPairAlignerThread* T = new ColorPairAlignerThread(seq_id_to_color_id,
                                                 writer,
                                                 1024*1024,
                                                 index,
                                                 *this,
                                                 cores);
            threads.push_back(T);
        }

        sbwt::SeqIO::Reader<> reader(fasta_file);
        run_dispatcher(threads, reader, 1024*1024);

        std::vector<std::int64_t> largest_color_ids;
        for (DispatcherConsumerCallback* t : threads) {
            ColorPairAlignerThread* cpat = static_cast<ColorPairAlignerThread*>(t);
            largest_color_ids.push_back(cpat->get_largest_color_id());
            delete t;
        }

        largest_color_id = *std::max_element(largest_color_ids.begin(), largest_color_ids.end());

        writer.flush();

        return outfile;
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

    void build_representation(const std::string& infile, const sdsl::bit_vector& cores, int64_t colorset_sampling_distance, int64_t ram_bytes, int64_t n_threads) {

        SBWT_backward_traversal_support backward_support(index_ptr);

        Buffered_ifstream<> in(infile, ios::binary);
        vector<char> buffer(16);

        vector<std::int64_t> node_set; // Reusable space
        vector<std::int64_t> colors_set; // Reusable space

        std::size_t set_id = 0;

        Sparse_Uint_Array_Builder builder(cores.size(), ram_bytes, n_threads);
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

            sets.emplace_back(colors_set);
            total_color_set_length += colors_set.size();

            for (const auto node : node_set) {
                builder.add(node, set_id);
                store_samples_in_unitig(cores, builder, backward_support, node, set_id, colorset_sampling_distance);
            }

            ++set_id;
        }

        node_id_to_color_set_id = builder.finish();

    }

    // Walks backward from from_node and marks every colorset_sampling_distance node on the way
    void store_samples_in_unitig(const sdsl::bit_vector& cores, Sparse_Uint_Array_Builder& builder, SBWT_backward_traversal_support& backward_support, int64_t from_node, int64_t colorset_id, int64_t colorset_sampling_distance){

        assert(cores[from_node] == 1);
        int64_t in_neighbors[4];
        int64_t indegree;
        backward_support.list_DBG_in_neighbors(from_node, in_neighbors, indegree);
        for(LL i = 0; i < indegree; i++){
            LL u = in_neighbors[i];
            LL counter = 0;
            while(cores[u] == 0){
                counter++;
                if(counter == colorset_sampling_distance){
                    builder.add(u, colorset_id);
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
    Coloring() {}

    Coloring(const std::vector<colorset_t>& sets,
             const Sparse_Uint_Array& node_id_to_color_set_id,
             const plain_matrix_sbwt_t& index) : sets(sets), node_id_to_color_set_id(node_id_to_color_set_id), index_ptr(&index){
    }

    std::size_t serialize(std::ostream& os) const {
        std::size_t bytes_written = 0;

        if(std::is_same<colorset_t, Bitmap_Or_Deltas_ColorSet>::value){
            string type_id = "bitmap-or-deltas-v1";
            bytes_written += sbwt::serialize_string(type_id, os);
        } else if(std::is_same<colorset_t, Roaring_Color_Set>::value){
            string type_id = "roaring-v0";
            bytes_written += sbwt::serialize_string(type_id, os);
        } else{
            throw std::runtime_error("Unsupported color set template");
        }

        std::size_t n_sets = sets.size();
        os.write(reinterpret_cast<char*>(&n_sets), sizeof(std::size_t));
        bytes_written += sizeof(std::size_t);

        for (std::size_t i = 0; i < n_sets; ++i) {
            bytes_written += sets[i].serialize(os);
        }

        bytes_written += node_id_to_color_set_id.serialize(os);

        os.write((char*)&largest_color_id, sizeof(largest_color_id));
        bytes_written += sizeof(largest_color_id);

        os.write((char*)&total_color_set_length, sizeof(total_color_set_length));
        bytes_written += sizeof(total_color_set_length);

        return bytes_written;
    }

    void load(std::ifstream& is, const plain_matrix_sbwt_t& index) {
        index_ptr = &index;

        string type_id = sbwt::load_string(is);

        // Check that the type id is correct for this class
        if(type_id == "bitmap-or-deltas-v1"){
            if(!std::is_same<colorset_t, Bitmap_Or_Deltas_ColorSet>::value){
                throw WrongTemplateParameterException();
            }
        } else if(type_id == "roaring-v0"){
            if(!std::is_same<colorset_t, Roaring_Color_Set>::value){
                throw WrongTemplateParameterException();
            }
        } else{
            throw std::runtime_error("Unknown color set type:" + type_id);
        }

        std::size_t n_sets = 0;
        is.read(reinterpret_cast<char*>(&n_sets), sizeof(std::size_t));

        sets.resize(n_sets);
        for (std::size_t i = 0; i < n_sets; ++i) {
            colorset_t cs;
            cs.load(is);
            sets[i] = cs;
        }

        node_id_to_color_set_id.load(is);

        is.read((char*)&largest_color_id, sizeof(largest_color_id));
        is.read((char*)&total_color_set_length, sizeof(total_color_set_length));
    }

    void load(const std::string& filename, const plain_matrix_sbwt_t& index) {
        throwing_ifstream in(filename, ios::binary);
        load(in.stream, index);
    }

    std::int64_t get_color_set_id(std::int64_t node) const {
        const auto& C_array = index_ptr->get_C_array();
        const auto& subset_struct= index_ptr->get_subset_rank_structure();

        while (!is_core_kmer(node)) {
            // While we don't have the color set id stored for the current node...

            // Follow an edge forward. The code below works only if we are at the
            // start of a suffix group. But this is guaranteed by the core k-mer marking
            // rules. If a suffix group is wider than 1, then all its elements are marked
            // as core because:
            //   - If there is at least one outgoing edge from the group, the nodes are marked
            //     by core k-mer rule (3) (see core_kmer_marker.hh)
            //   - If there are no outgoing edges from the group, the nodes are marked by
            //     core k-mer rule (2) (see core_kmer_marker.hh).
            if (subset_struct.A_bits[node] == 1) {
                node = C_array[0] + subset_struct.rank(node, 'A');
            } else if (subset_struct.C_bits[node] == 1) {
                node = C_array[1] + subset_struct.rank(node, 'C');
            } else if (subset_struct.G_bits[node] == 1) {
                node = C_array[2] + subset_struct.rank(node, 'G');
            } else if (subset_struct.T_bits[node] == 1) {
                node = C_array[3] + subset_struct.rank(node, 'T');
            } else {
                throw std::runtime_error("BUG: dead end in get_color_set_id");
            }
        }

        return node_id_to_color_set_id.get(node);
    }

    const colorset_t& get_color_set_of_node(std::int64_t node) const {
        std::int64_t color_set_id = get_color_set_id(node);
        return get_color_set_by_color_set_id(color_set_id);
    }

    // Yeah these function names are getting a bit verbose but I want to make it super clear
    // that the parameter is a color-set id and not a node id.
    const colorset_t& get_color_set_by_color_set_id(std::int64_t color_set_id) const {
        if (color_set_id == -1)
            throw std::runtime_error("BUG: Tried to access a color set with id " + to_string(color_set_id));
        return sets[color_set_id];
    }

    // Note! This function returns a new vector instead of a const-reference. Keep this
    // in mind if programming for performance. In that case, it's probably better to get the
    // color set using `get_color_set_of_node`, which returns a const-reference to a colorset_t object.
    std::vector<std::int64_t> get_color_set_of_node_as_vector(std::int64_t node) const {
        assert(node >= 0);
        assert(node < node_id_to_color_set_id.size());
        return get_color_set_of_node(node).get_colors_as_vector();
    }

    // See the comment on `get_color_set_of_node_as_vector`.
    std::vector<std::int64_t> get_color_set_as_vector_by_color_set_id(std::int64_t color_set_id) const {
        return get_color_set_by_color_set_id(color_set_id).get_colors_as_vector();
    }

    // If a node is a core k-mer, it has out-degree 1 and the color set of the out-neighbor is the
    // same as the color set of the node.
    bool is_core_kmer(std::int64_t node) const{
        return node_id_to_color_set_id.has_index(node);
    }

    void add_colors(const plain_matrix_sbwt_t& index,
                    const std::string fastafile,
                    const std::vector<std::int64_t>& colors_assignments,
                    const std::int64_t ram_bytes,
                    const std::int64_t n_threads,
                    int64_t colorset_sampling_distance) {

        index_ptr = &index;

        write_log("Marking core kmers", LogLevel::MAJOR);
        core_kmer_marker ckm;
        ckm.mark_core_kmers(fastafile, index);
        sdsl::bit_vector cores = ckm.core_kmer_marks;

        write_log("Getting node color pairs", LogLevel::MAJOR);
        const std::string node_color_pairs = get_node_color_pairs(index, fastafile, colors_assignments, cores, n_threads);

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
        build_representation(collected_nodes, cores, colorset_sampling_distance, ram_bytes, n_threads);
        get_temp_file_manager().delete_file(collected_nodes);

        write_log("Representation built", LogLevel::MAJOR);
    }

    int64_t largest_color() const{
        return largest_color_id;
    }

    int64_t number_of_distinct_color_sets() const{
        return sets.size();
    }

    int64_t sum_of_all_distinct_color_set_lengths() const{
        return total_color_set_length;
    }

    const std::vector<colorset_t>& get_all_distinct_color_sets() const{
        return sets;
    }

    // Returns map: component -> number of bytes
    map<string, int64_t> space_breakdown() const{
        map<string, int64_t> breakdown;
        int64_t color_set_total_size = 0;
        sbwt::SeqIO::NullStream ns;
        for(const colorset_t& cs : sets) color_set_total_size += cs.serialize(ns);
        breakdown["distinct-color-sets"] = color_set_total_size;
        

        for(auto [component, bytes] : node_id_to_color_set_id.space_breakdown()){
            breakdown["node-id-to-color-set-id-" + component] = bytes;
        }

        return breakdown;
    }

};


// Load whichever coloring data structure type is stored on disk
// The returned pointer must be eventually freed by the caller with delete
void load_coloring(string filename, const plain_matrix_sbwt_t& SBWT, 
std::variant<
Coloring<Bitmap_Or_Deltas_ColorSet>, 
Coloring<Roaring_Color_Set>>& coloring);