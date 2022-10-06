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

#include "roaring/roaring.hh"


class Color_Set {
    Roaring roaring;

public:
    Color_Set() {}

    Color_Set(Roaring r) : roaring(r) {}

    Color_Set(const vector<std::int64_t>& colors) {
        for (const auto x : colors)
            roaring.add(x);

        roaring.runOptimize();
    }

    Color_Set(const std::size_t n, const std::uint32_t* colors) {
        roaring.addMany(n, colors);

        roaring.runOptimize();
    }

    void add(const vector<std::int64_t>& colors) {
        for (const auto x : colors)
            roaring.add(x);

        roaring.runOptimize();
    }

    void add(const std::size_t n, const std::uint32_t* colors) {
        roaring.addMany(n, colors);

        roaring.runOptimize();
    }

    std::vector<std::uint32_t> get_colors_as_vector() const {
        std::vector<std::uint32_t> v(roaring.cardinality());
        roaring.toUint32Array(v.data());

        return v;
    }

    std::size_t size() const {
        return roaring.cardinality();
    }

    bool contains(const std::uint32_t n) const {
        return contains(n);
    }

    Color_Set intersection(const Color_Set& c) const {
        return Color_Set(roaring & c.roaring);
    }

    // union is a reserved word in C++ so this function is called do_union
    Color_Set do_union(const Color_Set& c) const {
        return Color_Set(roaring | c.roaring);
    }

    std::size_t serialize(std::ostream& os) const {
        std::size_t expected_size = roaring.getSizeInBytes();
        char* serialized_bytes = new char[expected_size];

        roaring.write(serialized_bytes);
        os.write(reinterpret_cast<char*>(&expected_size), sizeof(std::size_t));
        os.write(serialized_bytes, expected_size);
        delete[] serialized_bytes;

        return sizeof(std::size_t) + expected_size;
    }

    std::size_t load(std::ifstream& is) {
        std::size_t n;
        is.read(reinterpret_cast<char*>(&n), sizeof(std::size_t));

        char* serialized_bytes = new char[n];
        is.read(serialized_bytes, n);
        roaring = Roaring::readSafe(serialized_bytes, n);
        delete[] serialized_bytes;

        return sizeof(std::size_t) + n;
    }
};


class Coloring {
    std::vector<Color_Set> sets;
    
    Sparse_Uint_Array node_id_to_color_set_id;
    const plain_matrix_sbwt_t* index_ptr;

private:

    Coloring(const Coloring& temp_obj) = delete; // No copying
    Coloring& operator=(const Coloring& temp_obj) = delete;  // No copying

public:
    Coloring() {}

    Coloring(const std::vector<Color_Set>& sets,
             const sdsl::int_vector<>& color_set_pointers,
             const sdsl::bit_vector& cores,
             const plain_matrix_sbwt_t& index) : sets(sets) {
        index_ptr = &index;
    }

    std::size_t serialize(std::ostream& os) const {
        std::size_t bytes_written = 0;

        std::size_t n_sets = sets.size();
        os.write(reinterpret_cast<char*>(&n_sets), sizeof(std::size_t));
        bytes_written += sizeof(std::size_t);

        for (std::size_t i = 0; i < n_sets; ++i) {
            bytes_written += sets[i].serialize(os);
        }

        bytes_written += node_id_to_color_set_id.serialize(os);

        return bytes_written;
    }

    void load(std::ifstream& is, const plain_matrix_sbwt_t& index) {
        index_ptr = &index;

        std::size_t n_sets = 0;
        is.read(reinterpret_cast<char*>(&n_sets), sizeof(std::size_t));

        for (std::size_t i = 0; i < n_sets; ++i) {
            Color_Set cs;
            cs.load(is);
            sets.push_back(cs);
        }

        node_id_to_color_set_id.load(is);
    }

    void load(const std::string& filename, const plain_matrix_sbwt_t& index) {
        throwing_ifstream in(filename, ios::binary);
        load(in.stream, index);
    }

    std::int64_t get_color_set_id(std::int64_t node) const {
        const auto& C_array = index_ptr->get_C_array();
        const auto& subset_struct= index_ptr->get_subset_rank_structure();

        while (!node_id_to_color_set_id.has_index(node)) {
            // While we don't have the color set id stored for the current node...

            // Follow an edge forward
            if (subset_struct.A_bits[node] == 1) {
                node = C_array[0] + subset_struct.rank(node, 'A');
            } else if (subset_struct.C_bits[node] == 1) {
                node = C_array[1] + subset_struct.rank(node, 'C');
            } else if (subset_struct.G_bits[node] == 1) {
                node = C_array[2] + subset_struct.rank(node, 'G');
            } else if (subset_struct.T_bits[node] == 1) {
                node = C_array[3] + subset_struct.rank(node, 'T');
            } else {
                return -1;
            }
        }

        return node_id_to_color_set_id.get(node);
    }

    Color_Set get_color_set(std::int64_t node) const {
        std::int64_t color_set_id = get_color_set_id(node);

        if (color_set_id == -1) {
            Color_Set empty;
            return empty;
        }

        return sets[color_set_id];
    }

    std::vector<std::uint32_t> get_color_set_as_vector(std::int64_t node) const {
        assert(node >= 0);
        assert(node < node_id_to_color_set_id.size());
        return get_color_set(node).get_colors_as_vector();
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
        const std::string node_color_pairs = get_node_color_pairs(index, fastafile, colors_assignments, cores);

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

    std::string get_node_color_pairs(const plain_matrix_sbwt_t& index,
                                     const std::string& fasta_file,
                                     const std::vector<std::int64_t>& seq_id_to_color_id,
                                     const sdsl::bit_vector& cores) {
        const std::string outfile = get_temp_file_manager().create_filename();
        Buffered_ofstream<> out(outfile);
        const auto k = index.get_k();

        SeqIO::Reader<> reader(fasta_file);
        std::size_t seq_id = 0;
        std::size_t read_len = 0;

        while ((read_len = reader.get_next_read_to_buffer()) > 0) {

            if (read_len >= k) {
                const std::string seq(reader.read_buf, read_len);

                const auto res = index.streaming_search(seq);
                for (const auto node : res) {
                    if(node >= 0 && cores[node] == 1) {
                        write_big_endian_LL(out, node);
                        write_big_endian_LL(out, seq_id_to_color_id.at(seq_id));
                    }
                }
            }

            ++seq_id;

        }
        out.flush();

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
        
        string node_to_color_id_pairs_filename = get_temp_file_manager().create_filename();
        Buffered_ofstream<> node_to_color_id_pairs_out(node_to_color_id_pairs_filename, ios::binary);

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

            for (const auto node : node_set) {
                builder.add(node, set_id);
                store_samples_in_unitig(cores, builder, backward_support, node, set_id, colorset_sampling_distance, node_to_color_id_pairs_out);
            }

            ++set_id;
        }

        node_id_to_color_set_id = builder.finish();

    }

    // Walks backward from from_node and marks every colorset_sampling_distance node on the way
    void store_samples_in_unitig(const sdsl::bit_vector& cores, Sparse_Uint_Array_Builder& builder, SBWT_backward_traversal_support& backward_support, int64_t from_node, int64_t colorset_id, int64_t colorset_sampling_distance, Buffered_ofstream<>& out){

        int64_t n_new_marks = 0;
        
        assert(cores[from_node] == 1);
        int64_t in_neighbors[4];
        int64_t indegree;
        backward_support.list_in_neighbors(from_node, in_neighbors, indegree);
        for(LL i = 0; i < indegree; i++){
            LL u = in_neighbors[i];
            LL counter = 0;
            while(cores[u] == 0){
                counter++;
                if(counter == colorset_sampling_distance){
                    builder.add(u, colorset_id);
                    counter = 0;
                }
                backward_support.list_in_neighbors(u, in_neighbors, indegree);
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
    
};
