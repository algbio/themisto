#pragma once

#include "sbwt/SBWT.hh"
#include "sbwt/variants.hh"
#include "sdsl/select_support_mcl.hpp"

using namespace sbwt;

class SBWT_backward_traversal_support{

    private:

        // No copying
        SBWT_backward_traversal_support(throwing_ofstream const& other) = delete;
        SBWT_backward_traversal_support& operator=(throwing_ofstream const& other) = delete;

        const plain_matrix_sbwt_t* SBWT;  // Non-owning pointer
        sdsl::select_support_mcl<> select_A, select_C, select_G, select_T;
        sdsl::bit_vector dummy_marks;

        int64_t get_char_idx(char c) const{
            switch(c){
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T': return 3;
                default: return -1;
            }
        }

        char get_incoming_character(int64_t node) const{
            if(node < SBWT->get_C_array()[0]) return '$';
            if(node < SBWT->get_C_array()[1]) return 'A';
            if(node < SBWT->get_C_array()[2]) return 'C';
            if(node < SBWT->get_C_array()[3]) return 'G';
            return 'T';
        }

    public:

        SBWT_backward_traversal_support(const plain_matrix_sbwt_t* SBWT) : SBWT(SBWT){
            if(!SBWT->has_streaming_query_support())
                throw std::runtime_error("Bug: SBWT Streaming query support (=suffix group marks) required for backward traversal.");
            sdsl::util::init_support(select_A, &SBWT->get_subset_rank_structure().A_bits);
            sdsl::util::init_support(select_C, &SBWT->get_subset_rank_structure().C_bits);
            sdsl::util::init_support(select_G, &SBWT->get_subset_rank_structure().G_bits);
            sdsl::util::init_support(select_T, &SBWT->get_subset_rank_structure().T_bits);
            dummy_marks = SBWT->compute_dummy_node_marks();
        }

        const sdsl::bit_vector& get_dummy_marks() const{return dummy_marks;}

        // Up to 4 in-neighbors will be stored to the given array. The in-degree
        // will be stored to the other parameter. The in-degree means the in-degree in the
        // de Bruijn graph, not the SBWT graph (which has in-degree 1 everywhere except at the root).
        void list_DBG_in_neighbors(int64_t node, int64_t in_neighbors[4], int64_t& indegree) const{
            if(node == 0){ // Root
                indegree = 0; 
                return; 
            }

            node = backward_step(node);
            
            // This node and everything in the suffix group of this node are in-neighbors of the original node
            indegree = 0;
            if(dummy_marks[node] == 0) in_neighbors[indegree++] = node;
            node++;
            while(node < SBWT->number_of_subsets() && SBWT->get_streaming_support()[node] != 1){
                // Still in the same suffix group. This node can not be a dummy because there is
                // at most 1 dummy to each node.
                in_neighbors[indegree++] = node++;
            }          
        }

        // Goes one step backward in the SBWT. This is a well-defined operation because in the
        // SBWT each node except fot he root has exactly one incoming edge (even if the node has 
        // multiple in-neighbors in the DBG. If called on the root, we return back the root.
        int64_t backward_step(int64_t node) const{
            char c = get_incoming_character(node);
            if(c == '$') return node; // Root
            if(c == 'A') return select_A.select(node - SBWT->get_C_array()[0] + 1);
            if(c == 'C') return select_C.select(node - SBWT->get_C_array()[1] + 1);
            if(c == 'G') return select_G.select(node - SBWT->get_C_array()[2] + 1);
            if(c == 'T') return select_T.select(node - SBWT->get_C_array()[3] + 1);
            else throw std::runtime_error("This should never happen");
        }

        // Returns the incoming path label of length k to the node.
        // If the node is a dummy node and hence such path does not exist,
        // pads the label with dollars from the left.
        string get_node_label(int64_t node) const{
            int64_t k = SBWT->get_k();
            string label(k, '\0');
            for(int64_t i = 0; i < k; i++){
                char c = get_incoming_character(node);
                label[k-1-i] = c;
                if(c != '$'){
                    // Walk backward one step
                    if(c == 'A')
                        node = select_A.select(node - SBWT->get_C_array()[0] + 1);
                    if(c == 'C')
                        node = select_C.select(node - SBWT->get_C_array()[1] + 1);
                    if(c == 'G')
                        node = select_G.select(node - SBWT->get_C_array()[2] + 1);
                    if(c == 'T')
                        node = select_T.select(node - SBWT->get_C_array()[3] + 1);
                }
            }
            return label;
        }

};