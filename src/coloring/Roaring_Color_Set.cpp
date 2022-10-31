#include "coloring/Roaring_Color_Set.hh"

Roaring_Color_Set::Roaring_Color_Set(const Roaring_Color_Set::view_t& view){
    *this = *(view.ptr); // Make a copy
}