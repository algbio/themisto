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

        plain_matrix_sbwt_t* SBWT;  // Non-owning pointer
        sdsl::select_support_mcl<> select_A, select_C, select_G, select_T;

        char get_incoming_character(int64_t node){
            if(node < SBWT->get_C_array()[0]) return '$';
            if(node < SBWT->get_C_array()[1]) return 'A';
            if(node < SBWT->get_C_array()[2]) return 'C';
            if(node < SBWT->get_C_array()[3]) return 'G';
            return 'T';
        }

    public:

        SBWT_backward_traversal_support(plain_matrix_sbwt_t* SBWT) : SBWT(SBWT){
            if(!SBWT->has_streaming_query_support())
                throw std::runtime_error("Bug: SBWT Streaming query support (=suffix group marks) required for backward traversal.");
            sdsl::init_support(select_A, &SBWT->get_subset_rank_structure().A_bits);
            sdsl::init_support(select_C, &SBWT->get_subset_rank_structure().C_bits);
            sdsl::init_support(select_G, &SBWT->get_subset_rank_structure().G_bits);
            sdsl::init_support(select_T, &SBWT->get_subset_rank_structure().T_bits);
        }

        // Up to 4 in-neighbors will be stored to the given array. The in-degree
        // will be stored to the other parameter.
        void list_in_neighbors(int64_t node, int64_t in_neighbors[4], int64_t& indegree){
            char c = get_incoming_character(node);
            indegree = 0;
            if(c == '$') return; // Indegree 0

            // Walk backward one step
            if(c == 'A')
                node = select_A.select(node - SBWT->get_C_array()[0] + 1);
            if(c == 'C')
                node = select_C.select(node - SBWT->get_C_array()[1] + 1);
            if(c == 'G')
                node = select_G.select(node - SBWT->get_C_array()[2] + 1);
            if(c == 'T')
                node = select_T.select(node - SBWT->get_C_array()[3] + 1);
            
            // This node and everything in the suffix group of this node are in-neighbors of the original node
            indegree = 0;
            in_neighbors[indegree++] = node++;
            while(node < SBWT->number_of_subsets() && SBWT->get_streaming_support()[node] != 1){
                // Still in the same suffix group
                in_neighbors[indegree++] = node++;
            }
            
        }

};