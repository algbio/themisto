#pragma once

#include "globals.hh"
#include "sbwt/globals.hh"
#include "sbwt/variants.hh"
#include "coloring/Coloring.hh"
#include <vector>

using namespace std;

// Builds from existing index and serializes to disk
template<typename old_coloring_t, typename new_coloring_t> 
void build_from_index(plain_matrix_sbwt_t& dbg, const old_coloring_t& old_coloring, const string& to_index_dbg, const string& to_index_colors){

    // TODO: This makes a ton of unnecessary copies of things and has high peak RAM

    vector<typename new_coloring_t::colorset_type> new_colorsets;

    int64_t largest_color = 0;
    int64_t total_length = 0;

    for(int64_t i = 0; i < old_coloring.number_of_distinct_color_sets(); i++){
        vector<int64_t> set = old_coloring.get_color_set_as_vector_by_color_set_id(i);

        for(int64_t x : set) largest_color = max(largest_color, x);
        total_length += set.size();

        new_colorsets.push_back(set); // Re-encodes as new_coloring_t::color_set_t
    }

    typename new_coloring_t::colorset_storage_type new_storage(new_colorsets);
    new_coloring_t new_coloring(new_storage, 
                                old_coloring.get_node_id_to_colorset_id_structure(),
                                dbg, largest_color, total_length);

    write_log("Serializing to " + to_index_dbg + " and " + to_index_colors, LogLevel::MAJOR);

    throwing_ofstream colors_out(to_index_colors);
    new_coloring.serialize(colors_out.stream);
    
    throwing_ofstream dbg_out(to_index_dbg);
    dbg.serialize(dbg_out.stream);
}

void transform_existing_index(const string& from_index_dbg, const string& from_index_coloring, const string& to_index_dbg, const string& to_index_coloring, const string& new_index_color_set_type){

    write_log("Building new structure of type " + new_index_color_set_type, LogLevel::MAJOR);

    // Transform existing index to new format

    sbwt::write_log("Loading de Bruijn Graph", sbwt::LogLevel::MAJOR);
    std::unique_ptr<sbwt::plain_matrix_sbwt_t> dbg_ptr = std::make_unique<sbwt::plain_matrix_sbwt_t>();
    dbg_ptr->load(from_index_dbg);

    sbwt::write_log("Loading coloring", sbwt::LogLevel::MAJOR);
    std::variant<Coloring<SDSL_Variant_Color_Set>, Coloring<Roaring_Color_Set>> old_coloring;
    load_coloring(from_index_coloring, *dbg_ptr, old_coloring);

    if(std::holds_alternative<Coloring<SDSL_Variant_Color_Set>>(old_coloring))
        write_log("sdsl coloring structure loaded", LogLevel::MAJOR);
    if(std::holds_alternative<Coloring<Roaring_Color_Set>>(old_coloring))
        write_log("roaring coloring structure loaded", LogLevel::MAJOR);

    auto visitor = [&](auto& old){
        if(new_index_color_set_type == "sdsl-hybrid"){
            build_from_index<decltype(old), Coloring<SDSL_Variant_Color_Set>>(*dbg_ptr, old, to_index_dbg, to_index_coloring);
        } else if(new_index_color_set_type == "roaring"){
            build_from_index<decltype(old), Coloring<Roaring_Color_Set>>(*dbg_ptr, old,  to_index_dbg, to_index_coloring);
        } else{
            throw std::runtime_error("Unkown coloring structure type: " + new_index_color_set_type);
        }
    };

    std::visit(visitor, old_coloring);

    return;

}