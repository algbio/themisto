#include "KallistoLite.hh"
#include "input_reading.hh"
#include <string>
#include <cstring>

using namespace std;

LL count_ones(sdsl::bit_vector& bv){
    LL ones = 0;
    for(LL i = 0; i < bv.size(); i++) ones += bv[i];
    return ones;
}

LL count_zeros(sdsl::bit_vector& bv){
    return bv.size() - count_ones(bv);
}

int main2(int argc, char** argv){

    if(argc == 1){
        cerr << "Usage: ./program index_dir" << endl;
        cerr << "Writes to stdout" << endl;
        exit(1);
    }

    string index_dir = argv[1];
    check_dir_exists(index_dir);

    KallistoLite kl;
    kl.load_boss(index_dir + "/boss-");
    kl.load_colors(index_dir + "/coloring-");

    LL n_nodes;
    LL n_edges;
    LL n_branching_nodes;
    LL n_colors;
    LL n_nonempty_colorsets;
    LL n_nonredundant_colorsets;
    LL n_nonempty_and_nonredundant_colorsets;

    n_nodes = kl.boss.get_number_of_nodes();
    n_edges = kl.boss.get_number_of_edges();
    n_branching_nodes = 0;
    for(LL i = 0; i < n_nodes; i++){
        n_branching_nodes += kl.boss.is_branching(i);
    }
    n_colors = kl.coloring.n_colors;
    n_nonempty_colorsets = count_ones(kl.coloring.nonempty);
    n_nonredundant_colorsets = count_zeros(kl.coloring.redundancy_marks);
    n_nonempty_and_nonredundant_colorsets = count_ones(kl.coloring.nonempty_and_nonredundant);

    cout << "number of nodes: " << n_nodes << endl;
    cout << "number of edges: " << n_edges << endl;
    cout << "branching nodes: " << n_branching_nodes << " (" << (double)n_branching_nodes / n_nodes << ")" << endl;
    cout << "number of colors: " << n_colors << endl;
    cout << "number of nonempty colorsets: " << n_nonempty_colorsets << " (" << 
             (double)n_nonempty_colorsets / n_nodes << ")" << endl;
    cout << "number of nonredundant colorsets: " << n_nonredundant_colorsets << " (" <<
             (double)n_nonredundant_colorsets / n_nodes << ")" << endl;
    cout << "nonempty and nonredundant: " << n_nonempty_and_nonredundant_colorsets << " (" <<
          (double) n_nonempty_and_nonredundant_colorsets / n_nodes << ")" << endl;

    if(n_colors > (1LL << 31)) return 0; // Does not fit in a 32-bit integer

    LL n_distinct_color_sets;
    set<vector<int32_t> > colorsets;
    for(LL i = 0; i < n_nodes; i++){
        set<LL> S = kl.coloring.get_colorset(i, kl.boss);
        vector<int32_t> colorset(S.begin(), S.end());
        colorsets.insert(colorset);
    }

    n_distinct_color_sets = colorsets.size();

    cout << "number of distinct colorsets " << n_distinct_color_sets << " (" <<
            (double)n_distinct_color_sets / n_nodes << ")" << endl;


    // Number of nodes
    // Number of edges
    // Number of redundant nodes
    // Number of branching nodes
    // Number of colors
    // Number of distinct color sets

}

int main(int argc, char** argv){
    try{
        return main2(argc, argv);
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    }
}