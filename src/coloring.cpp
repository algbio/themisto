#include "Coloring.hh"

void load_coloring(string filename, const plain_matrix_sbwt_t& SBWT,
std::variant<
Coloring<SDSL_Variant_Color_Set>,
Coloring<Roaring_Color_Set>>& coloring){

    Coloring<SDSL_Variant_Color_Set> coloring1;
    Coloring<Roaring_Color_Set> coloring2;

    try{
        throwing_ifstream colors_in(filename, ios::binary);
        coloring = coloring1;
        std::get<Coloring<SDSL_Variant_Color_Set>>(coloring).load(colors_in.stream, SBWT);
        return; // No exception thrown
    } catch(Coloring<SDSL_Variant_Color_Set>::WrongTemplateParameterException& e){
        // Was not this one
    }

    try{
        throwing_ifstream colors_in(filename, ios::binary);
        coloring = coloring2;
        std::get<Coloring<Roaring_Color_Set>>(coloring).load(colors_in.stream, SBWT);
        return; // No exception thrown
    } catch(Coloring<Roaring_Color_Set>::WrongTemplateParameterException& e){
        // Was not this one
    }

    throw std::runtime_error("Error: could not load color structure.");
}
