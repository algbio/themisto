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
#include "Fixed_Width_Int_Color_Set.hh"
#include <variant>
#include "bit_magic_color_set.hh"

template <typename T>
concept Color_Set_Interface = requires(T& t, std::ostream& os, std::istream& is){
    { t.empty() } -> std::same_as<bool>; // This should take constant time
    { t.size() } -> std::same_as<int64_t>; // This may take linear time
    { t.size_in_bits() } -> std::same_as<int64_t>;
    { t.contains(int64_t()) } -> std::same_as<bool>; // This may take linear time
    { t.intersection(t) } -> std::same_as<T>;
    { t.do_union(t) } -> std::same_as<T>;
    { t.serialize(os) } -> std::same_as<int64_t>; // Returns the number of bytes written
    { t.load(is) } -> std::same_as<void>;
    { t.get_colors_as_vector() } -> std::same_as<std::vector<int64_t>>;
    requires std::constructible_from<T, vector<int64_t>>;
    requires std::default_initializable<T>;
};

// Takes as parameter a class that encodes a single color set
template<typename colorset_t = Bitmap_Or_Deltas_ColorSet> requires Color_Set_Interface<colorset_t>
class Coloring {

private:    

    std::vector<colorset_t> sets;
    Sparse_Uint_Array node_id_to_color_set_id;
    const plain_matrix_sbwt_t* index_ptr;
    int64_t largest_color_id = 0;
    int64_t total_color_set_length = 0;

public:
    class WrongTemplateParameterException : public std::exception{
        const char * what() const noexcept override{
            return "Template type id in a serialized Coloring structure does not match the class template parameter.";
        }
    };

    typedef colorset_t colorset_type;

public:
    Coloring() {}

    Coloring(const std::vector<colorset_t>& sets,
             const Sparse_Uint_Array& node_id_to_color_set_id,
             const plain_matrix_sbwt_t& index) : sets(sets), node_id_to_color_set_id(node_id_to_color_set_id), index_ptr(&index){
    }

    std::size_t serialize(std::ostream& os) const {
        std::size_t bytes_written = 0;

        if(std::is_same<colorset_t, Fixed_Width_Int_Color_Set>::value){
            string type_id = "sdsl-fixed-v0";
            bytes_written += sbwt::serialize_string(type_id, os);
        } else if(std::is_same<colorset_t, Bitmap_Or_Deltas_ColorSet>::value){
            string type_id = "sdsl-hybrid-v1";
            bytes_written += sbwt::serialize_string(type_id, os);
        } else if(std::is_same<colorset_t, Roaring_Color_Set>::value){
            string type_id = "roaring-v0";
            bytes_written += sbwt::serialize_string(type_id, os);
        } else if(std::is_same<colorset_t, Bit_Magic_Color_Set>::value){
            string type_id = "bitmagic-v0";
            bytes_written += sbwt::serialize_string(type_id, os);
        }  else{
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
        if(type_id == "sdsl-fixed-v0"){
            if(!std::is_same<colorset_t, Fixed_Width_Int_Color_Set>::value){
                throw WrongTemplateParameterException();
            }
        }  else if(type_id == "sdsl-hybrid-v1"){
            if(!std::is_same<colorset_t, Bitmap_Or_Deltas_ColorSet>::value){
                throw WrongTemplateParameterException();
            }
        } else if(type_id == "roaring-v0"){
            if(!std::is_same<colorset_t, Roaring_Color_Set>::value){
                throw WrongTemplateParameterException();
            }
        } else if(type_id == "bitmagic-v0"){
            if(!std::is_same<colorset_t, Bit_Magic_Color_Set>::value){
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

    template<class Tprime> requires Color_Set_Interface<Tprime>
    friend class Coloring_Builder;
};

// Load whichever coloring data structure type is stored on disk
// The returned pointer must be eventually freed by the caller with delete
void load_coloring(string filename, const plain_matrix_sbwt_t& SBWT,
std::variant<
Coloring<Bitmap_Or_Deltas_ColorSet>,
Coloring<Roaring_Color_Set>,
Coloring<Fixed_Width_Int_Color_Set>,
Coloring<Bit_Magic_Color_Set>>& coloring);

