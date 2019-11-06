#pragma once

#include <set>
#include <unordered_map>
#include <cassert>
#include "globals.hh"


// Takes in a file with a list of names, one on each line.
// Assigns every distinct name a distinct id, starting from
// zero such that the new ids are handed out in the order of 
// appearance of the names in the file

class NameMapping{

public:

    unordered_map<LL, string> id_to_name;
    unordered_map<string, LL> name_to_id;

    NameMapping() {}
    NameMapping(string filename){
        build_from_namefile(filename);
    }

    // Color file: one line for each sequence: the name of the color (a string)
    void build_from_namefile(string filename){
        vector<string> naming = read_namefile(filename);
        build_from_name_vector(naming);
    }

    // Input: a vector where line i has the color of sequence i
    void build_from_name_vector(vector<string> naming){
        id_to_name.clear();
        name_to_id.clear();

        // Build the mapping
        name_to_id = assign_ids(naming);

        // Build the inverse mapping
        for(auto keyvalue : name_to_id){
            id_to_name[keyvalue.second] = keyvalue.first;
        }
    }

    LL get_number_of_names(){
        return name_to_id.size();
    }

    vector<LL> names_to_ids(vector<string> names){
        vector<LL> res;
        for(string name : names) res.push_back(name_to_id[name]);
        return res;
    }

    vector<string> ids_to_names(vector<LL> ids){
        vector<string> res;
        for(LL id : ids) res.push_back(id_to_name[id]);
        return res;
    }

    vector<string> ids_to_names(set<LL> ids){
        vector<string> res;
        for(LL id : ids) res.push_back(id_to_name[id]);
        return res;
    }

    void save_to_disk(string path_prefix){
        check_writable(path_prefix + "id_to_name");
        throwing_ofstream id_to_name_out(path_prefix + "id_to_name");
        for(auto keyvalue : id_to_name){
            id_to_name_out << keyvalue.first << " " << keyvalue.second << "\n";
        }
    }

    void load_from_disk(string path_prefix){ 
        id_to_name.clear();
        name_to_id.clear();

        check_readable(path_prefix + "id_to_name");    
        throwing_ifstream input(path_prefix + "id_to_name");
        string line;
        while(input.getline(line)){
            vector<string> tokens = parse_tokens<string>(line);
            assert(tokens.size() == 2);
            LL id = stoll(tokens[0]);
            string name = tokens[1];
            id_to_name[id] = name;
            name_to_id[name] = id;
        }
    }

private:

    vector<string> read_namefile(string namefile){
        throwing_ifstream namestream(namefile);

        vector<string> naming;
        string name;
        while(namestream.getline(name)){
            naming.push_back(name);
        }

        return naming;
    }

    unordered_map<string, LL> assign_ids(vector<string>& naming){
        // Assign ids in the order the colors first appear in the coloring
        LL id = 0;
        unordered_map<string, LL> mapping;
        for(string name : naming){
            if(mapping.find(name) == mapping.end()){
                mapping[name] = id;
                id++;
            }
        }
        return mapping;
    }


};

