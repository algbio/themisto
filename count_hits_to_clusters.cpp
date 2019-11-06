#include <string>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <algorithm>
#include "globals.hh"
#include "NameMapping.hh"

using namespace std;

/*

Assumes:

  - The color ids in the alignment results are the ref seq ids directly

*/

// Returns mapping ref id -> cluster id
// ColorMapping defines cluster name <-> cluster id
vector<LL> read_clusterfile(string clusterfile){
    NameMapping nm(clusterfile); // cluster id <-> cluster name
    throwing_ifstream stream(clusterfile);
    vector<LL> ref_id_to_cluster_id; // We are building this
    string cluster_name;
    while(stream.getline(cluster_name)){
        LL cluster_id = nm.name_to_id[cluster_name];
        ref_id_to_cluster_id.push_back(cluster_id);
    }

    return ref_id_to_cluster_id;
}

// Write pairs (cluster id, cluster name)
void write_mapping(string clusterfile, string outfile){
    throwing_ofstream output(outfile);
    NameMapping nm(clusterfile); // cluster id <-> cluster name
    for(LL i = 0; i < nm.get_number_of_names(); i++){
        output << i << " " << nm.id_to_name[i] << "\n";
    }
}

void write_hit_counts(string alignmentfile, string clusterfile, string outfile){

    // Need: color id -> cluster id
    // For the user, also need to write cluster id -> cluster name

    // Step 1: give distinct clusters ids
    // Step 2: iterate alignments

    // Assign clusters ids
    
    vector<LL> ref_id_to_cluster_id = read_clusterfile(clusterfile);

    throwing_ifstream input(alignmentfile);
    throwing_ofstream output(outfile);

    // For each read
    string line;
    while(input.getline(line)){
        unordered_map<LL, LL> hits; // cluster id -> number of hits
        vector<LL> numbers = parse_tokens<LL>(line);
        LL read_id = numbers[0];
        for(LL i = 1; i < numbers.size(); i++){
            // numbers[i] is a color id. Want a cluster id
            LL ref_id = numbers[i];
            LL cluster_id = ref_id_to_cluster_id[ref_id];
            hits[cluster_id]++;
        }

        output << read_id << " ";
        for(auto& keyvalue : hits){
            LL cluster_id = keyvalue.first;
            LL count = keyvalue.second;
            output << cluster_id << " " << count << " ";
        }
        output << "\n";
    }
}



int main2(int argc, char** argv){

    if(argc != 5){
        cerr << "Usage: ./program alignmentfile clusterfile counts-out mapping-out" << endl;
        exit(1);
    }

    string alignmentfile = argv[1];
    string clusterfile = argv[2];
    string counts_outfile = argv[3];
    string mapping_outfile = argv[4];

    check_readable(alignmentfile);
    check_readable(clusterfile);
    check_writable(counts_outfile);
    check_writable(mapping_outfile);

    write_mapping(clusterfile, mapping_outfile);
    write_hit_counts(alignmentfile, clusterfile, counts_outfile);


}

int main(int argc, char** argv){
    try{
        return main2(argc, argv);
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    }
}
