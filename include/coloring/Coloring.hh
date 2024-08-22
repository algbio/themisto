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
#include "SeqIO/buffered_streams.hh"
#include "sbwt/EM_sort/bit_level_stuff.hh"
#include "sbwt/EM_sort/EM_sort.hh"
#include "sbwt/globals.hh"
#include "SeqIO/SeqIO.hh"

#include <vector>

#include <cstdint>

#include "WorkDispatcher.hh"
#include "hybrid_color_set.hh"
#include "Roaring_Color_Set.hh"
#include "Color_Set_Storage.hh"
#include "Color_Set.hh"
#include "Color_Set_Interface.hh"
#include <variant>

// Takes as parameter a class that encodes a single color set
template<typename colorset_t = SDSL_Variant_Color_Set> 
requires Color_Set_Interface<colorset_t>
class Coloring {

public:

typedef colorset_t colorset_type;
typedef colorset_t::view_t colorset_view_type;
typedef Color_Set_Storage<colorset_t> colorset_storage_type;

class WrongTemplateParameterException : public std::exception{
    const char * what() const noexcept override{
        return "Template type id in a serialized Coloring structure does not match the class template parameter.";
    }
};

private:    

    colorset_storage_type sets;
    Sparse_Uint_Array node_id_to_color_set_id;
    const plain_matrix_sbwt_t* index_ptr;
    int64_t largest_color_id = 0;
    int64_t total_color_set_length = 0;

public:

    Coloring() {}

    Coloring(const colorset_storage_type& sets,
             const Sparse_Uint_Array& node_id_to_color_set_id,
             const plain_matrix_sbwt_t& index,
             const int64_t largest_id,
             const int64_t total_color_set_length) 
             : sets(sets), node_id_to_color_set_id(node_id_to_color_set_id), index_ptr(&index), largest_color_id(largest_id), total_color_set_length(total_color_set_length){
    }

    std::size_t serialize(std::ostream& os) const {
        std::size_t bytes_written = 0;

        if(std::is_same<colorset_t, SDSL_Variant_Color_Set>::value){
            string type_id = "sdsl-hybrid-v4";
            bytes_written += sbwt::serialize_string(type_id, os);
        } else if(std::is_same<colorset_t, Roaring_Color_Set>::value){
            string type_id = "roaring-v0";
            bytes_written += sbwt::serialize_string(type_id, os);
        } else{
            throw std::runtime_error("Unsupported color set template");
        }

        bytes_written += sets.serialize(os);
        bytes_written += node_id_to_color_set_id.serialize(os);

        os.write((char*)&largest_color_id, sizeof(largest_color_id));
        bytes_written += sizeof(largest_color_id);

        os.write((char*)&total_color_set_length, sizeof(total_color_set_length));
        bytes_written += sizeof(total_color_set_length);

        return bytes_written;
    }

    int64_t serialize(const string& filename) const{
        throwing_ofstream out(filename, ios::binary);
        return serialize(out.stream);
    }


    void load(std::ifstream& is, const plain_matrix_sbwt_t& index) {
        index_ptr = &index;

        string type_id = sbwt::load_string(is);

        // Check that the type id is correct for this class
        if(type_id == "sdsl-hybrid-v4"){
            if(!std::is_same<colorset_t, SDSL_Variant_Color_Set>::value){
                throw WrongTemplateParameterException();
            }
        } else if(type_id == "roaring-v0"){
            if(!std::is_same<colorset_t, Roaring_Color_Set>::value){
                throw WrongTemplateParameterException();
            }
        } else{
            throw std::runtime_error("Unknown color set type:" + type_id);
        }

        sets.load(is);
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

    colorset_view_type get_color_set_of_node(std::int64_t node) const {
        std::int64_t color_set_id = get_color_set_id(node);
        return get_color_set_by_color_set_id(color_set_id);
    }

    // Yeah these function names are getting a bit verbose but I want to make it super clear
    // that the parameter is a color-set id and not a node id.
    colorset_view_type get_color_set_by_color_set_id(std::int64_t color_set_id) const {
        if (color_set_id == -1)
            throw std::runtime_error("BUG: Tried to access a color set with id " + to_string(color_set_id));
        return sets.get_color_set_by_id(color_set_id);
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
        return sets.number_of_sets_stored();
    }

    int64_t sum_of_all_distinct_color_set_lengths() const{
        return total_color_set_length;
    }

    const Sparse_Uint_Array& get_node_id_to_colorset_id_structure() const{
        return node_id_to_color_set_id;
    }

    const std::vector<colorset_view_type> get_all_distinct_color_sets() const{
        return sets.get_all_sets();
    }

    // Returns map: component -> number of bytes
    map<string, int64_t> space_breakdown() const{
        map<string, int64_t> breakdown;

        for(auto [component, bytes] : sets.space_breakdown()){
            breakdown["color-set-storage-" + component] = bytes;
        }

        for(auto [component, bytes] : node_id_to_color_set_id.space_breakdown()){
            breakdown["node-id-to-color-set-id-" + component] = bytes;
        }

        return breakdown;
    }

    // Increases the index size, but makes queries faster
    void add_all_node_id_to_color_set_id_pointers(const plain_matrix_sbwt_t& index, SBWT_backward_traversal_support& sbwt_bws, int64_t n_threads) {

        // Data structure for the new "sparse" array of values
        uint64_t max_value = node_id_to_color_set_id.get_max_value();
        sdsl::bit_vector marks(index.number_of_subsets(), 1); // Everything is marked, even dummies, but that is ok
        sdsl::int_vector<> values(index.number_of_subsets(), 0, std::bit_width(max_value));

        // The code below is parallel, so the for-loop is in batches so that the critical
        // section is run less often
        int64_t batch_size = 10000;

        #pragma omp parallel for num_threads (n_threads)
        for(int64_t b = 0; b < index.number_of_subsets(); b += batch_size){
            int64_t batch_end = min(b + batch_size, index.number_of_subsets()); // One past the end
            vector<pair<int64_t,int64_t>> updates; // Batched updates

            for(int64_t v = b; v < batch_end; v++){

                int64_t in_neighbors[4];
                int64_t indegree;

                if(this->node_id_to_color_set_id.has_index(v)){
                    
                    int64_t value = this->node_id_to_color_set_id.get(v);
                    updates.push_back({v,value});

                    sbwt_bws.list_DBG_in_neighbors(v, in_neighbors, indegree);
                    for(int64_t i = 0; i < indegree; i++){
                        int64_t u = in_neighbors[i];
                        while(!this->node_id_to_color_set_id.has_index(u)){
                            updates.push_back({u,value});

                            sbwt_bws.list_DBG_in_neighbors(u, in_neighbors, indegree);
                            if(indegree == 0) break; // Root node
                            if(indegree >= 2) break; // Predecessors are already marked
                            u = in_neighbors[0]; // The only in-neighbor
                        }
                    }
                }
            }

            // Critical section: apply the updates
            #pragma omp critical
            {
                for(auto [v, value] : updates){
                    values[v] = value;
                }
            }
        }
        
        this->node_id_to_color_set_id = Sparse_Uint_Array(marks, values, max_value);

    }

    const plain_matrix_sbwt_t& get_sbwt() const {
        return *(this->index_ptr);
    }

    template<typename T1, typename T2> requires Color_Set_Interface<T1>
    friend class Coloring_Builder;

    template<typename T1> requires Color_Set_Interface<T1>
    friend class Coloring_Builder_From_GGCAT;
};

// Load whichever coloring data structure type is stored on disk
// The returned pointer must be eventually freed by the caller with delete
void load_coloring(string filename, const plain_matrix_sbwt_t& SBWT,
std::variant<
Coloring<SDSL_Variant_Color_Set>,
Coloring<Roaring_Color_Set>>& coloring);

