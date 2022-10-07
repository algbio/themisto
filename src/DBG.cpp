#include "DBG.hh"

DBG::all_nodes_generator DBG::all_nodes() const{
    return all_nodes_generator(SBWT, &(backward_support.get_dummy_marks()));
}

DBG::outedge_generator DBG::outedges(Node v) const{
    return outedge_generator(v, SBWT);
}

DBG::inedge_generator DBG::inedges(Node v) const{
    return inedge_generator(v, SBWT, &backward_support);
}
