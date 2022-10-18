#pragma once

#include <bitset>
#include <string>

#include <sdsl/bit_vectors.hpp>

#include "sbwt/SBWT.hh"
#include "SeqIO.hh"
#include "sbwt/variants.hh"
#include "globals.hh"

using namespace sbwt;

/*
    Marks nodes v that fall into one of the following four cases:

    (1) v has an outgoing edge into a node that is the first k-mer of a sequence
    (2) v is the last k-mer of a sequence
    (3) v has an outgoing edge into a node with in-degree at least 2
    (4) v has outdegree at least 2

    The idea of these cases is that a node u has the same color set as the nearest
    marked node v that can be reached by walking forward from u in the de Bruijn graph.
    In other words, the color set of a non-marked node u is the same as the color set of
    its successor in the graph (the successor exists due to case 2 and is unique due to case 4).
*/

template<typename sequence_reader_t = SeqIO::Reader<>>
class core_kmer_marker {
    static constexpr char int_to_dna[] = { 'A', 'C', 'G', 'T' };
public:
    sdsl::bit_vector core_kmer_marks;

    core_kmer_marker() {}

    core_kmer_marker(const std::size_t sz) : core_kmer_marks(sz, 0) {}

    std::size_t mark_core_kmers(sequence_reader_t& reader, const plain_matrix_sbwt_t& index, bool reverse_complements) {
        std::size_t total_core_count = 0;
        sdsl::util::assign(core_kmer_marks, sdsl::bit_vector(index.number_of_subsets(), 0));

        write_log("Handling cases one and two", LogLevel::MAJOR);
        total_core_count += handle_case_one_and_two(reader, index, reverse_complements);
        write_log("Handling case three", LogLevel::MAJOR);
        total_core_count += handle_case_three(index);
        write_log("Handling case four", LogLevel::MAJOR);
        total_core_count += handle_case_four(index);

        return total_core_count;
    }

    inline std::size_t handle_case_one_and_two(sequence_reader_t& reader, const plain_matrix_sbwt_t& index, bool reverse_complements){
        const std::size_t n = index.number_of_subsets();
        const auto k = index.get_k();
        const auto& rank_structure = index.get_subset_rank_structure();
        const auto& C_array = index.get_C_array();

        std::size_t cores = 0;

        sdsl::bit_vector first_kmer_marks(n, 0);
        std::size_t read_len = 0;
        while ((read_len = reader.get_next_read_to_buffer()) > 0) {

            vector<string> parts = split_at_non_ACGT(reader.read_buf, read_len);
            if(reverse_complements){
                int64_t n_parts = parts.size();
                for(int64_t i = 0; i < n_parts; i++){
                    parts.push_back(get_rc(parts[i]));
                }
            }

            for(const string& part : parts){
                if (part.size() >= k) {
                    // End of a sequence
                    const std::size_t last_kmer_idx = index.search(part.substr(part.size() - k));
                    core_kmer_marks[last_kmer_idx] = 1;
                    cores += core_kmer_marks[last_kmer_idx] == 1 ? 0 : 1;

                    // Beginning of a sequence
                    const std::size_t first_kmer_idx = index.search(part.substr(0, k));
                    first_kmer_marks[first_kmer_idx] = 1;
                }
            }
        }

        // Mark nodes preceding beginning of a sequence
        for (std::size_t i = 1; i < n; ++i) {
            const auto& edges = column_edges(index, i);
            if (edges.count() > 0) {
                for (std::size_t j = 0; j < 4; ++j) {
                    if (edges.test(j)) {
                        const std::size_t destination = C_array[j] + rank_structure.rank(i, int_to_dna[j]);

                        if (first_kmer_marks[destination] == 1) {
                            cores += core_kmer_marks[i] == 1 ? 0 : 1;
                            core_kmer_marks[i] = 1;
                        }
                    }
                }
            }
        }

        return cores;
    }

    inline std::size_t handle_case_three(const plain_matrix_sbwt_t& index) {
        const std::size_t n = index.number_of_subsets();
        const auto& suffix_group_marks = index.get_streaming_support();

        std::size_t cores = 0;

        for (std::size_t i = 1; i < n; ++i) {
            if (suffix_group_marks[i] == 1) {
                const std::size_t group_width = suffix_group_width(suffix_group_marks, i);

                if (group_width > 1) {
                    for (std::size_t j = 0; j < group_width; ++j) {
                        cores += core_kmer_marks[i + j] == 1 ? 0 : 1;
                        core_kmer_marks[i + j] = 1;
                    }
                }
            }
        }

        return cores;
    }

    inline std::size_t suffix_group_width(const sdsl::bit_vector& marks, const std::size_t idx) const {
        std::size_t width = 1;

        while (idx + width < marks.size() && marks[idx + width] != 1)
            ++width;

        return width;
    }

    inline std::size_t handle_case_four(const plain_matrix_sbwt_t& index) {
        const std::size_t n = index.number_of_subsets();

        std::size_t cores = 0;

        for (std::size_t i = 1; i < n; ++i) {
            const auto& edges = column_edges(index, i);

            // If node has a branch, mark it
            if (edges.count() > 1) {
                cores += core_kmer_marks[i] == 1 ? 0 : 1;
                core_kmer_marks[i] = 1;
            }
        }

        return cores;
    }

    inline std::bitset<4> column_edges(const plain_matrix_sbwt_t& index, const std::size_t idx) const {
        const auto& subset_struct = index.get_subset_rank_structure();
        std::bitset<4> edges{0b0000};

        if (subset_struct.A_bits[idx] == 1)
            edges.set(0, true);
        if (subset_struct.C_bits[idx] == 1)
            edges.set(1, true);
        if (subset_struct.G_bits[idx] == 1)
            edges.set(2, true);
        if (subset_struct.T_bits[idx] == 1)
            edges.set(3, true);

        return edges;
    }

    vector<string> split_at_non_ACGT(const char* S, int64_t S_size){
        vector<string> parts;
        int64_t start = 0;
        for(int64_t end = 0; end <= S_size; end++){ // Exclusive end point
            if(end == S_size || (S[end] != 'A' && S[end] != 'C' && S[end] != 'G' && S[end] != 'T')){
                if(end > start) {
                    parts.push_back(string(S + start, S + end));
                }
                start = end + 1;
            }
        }
        return parts;
    }

};
