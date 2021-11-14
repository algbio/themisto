#pragma once

#include <string>
#include <vector>
#include <fstream>
#include "globals.hh"
#include "EM_sort.hh"
#include "EM_algorithms.hh"

string EM_sort_big_endian_LL_pairs(string infile, LL ram_bytes, LL key, LL n_threads){
    assert(key == 0 || key == 1);
    
    auto cmp = [&](const char* A, const char* B) -> bool{
        LL x_1,y_1,x_2,y_2; // compare (x_1, y_1) and (x_2,y_2)
        x_1 = parse_big_endian_LL(A + 0);
        y_1 = parse_big_endian_LL(A + 8);
        x_2 = parse_big_endian_LL(B + 0);
        y_2 = parse_big_endian_LL(B + 8);
        if(key == 0){
            return make_pair(x_1, y_1) < make_pair(x_2, y_2);
        } else{
            return make_pair(y_1, x_1) < make_pair(y_2, x_2);
        }
    };

    string outfile = get_temp_file_manager().create_filename();
    EM_sort_constant_binary(infile, outfile, cmp, ram_bytes, 4, 8+8, n_threads);
    return outfile;
}

// Input: file with one ';'-separated pair of strings (x,y) on each line
//        cmp_1 compares x values and cmp_y compares y values. 
//        The comparison functions should return < 0, 0 or > 0 like strcmp.
// Output: a file with the lines sorted by x if key = 0, or by y if key == 1
//         if the key values are equal, the other elements are compared
string EM_line_sort_pairs(string infile, int (*cmp_1) (const char* x_1, const char* x_2), 
                                         int (*cmp_2) (const char* x_1, const char* x_2), 
                                         LL key, LL ram_bytes, LL n_threads){
    assert(key == 0 || key == 1);

    auto cmp = [&](const char* A, const char* B) -> bool{
        vector<string> splitted_A = split(A,';');
        vector<string> splitted_B = split(B,';');
        LL b = cmp_1(splitted_A[key].c_str(), splitted_B[key].c_str());
        if(b == 0) b = cmp_2(splitted_A[1-key].c_str(), splitted_B[1-key].c_str());
        return b < 0;
    };
    
    string outfile = get_temp_file_manager().create_filename();
    EM_sort(infile, outfile, cmp, ram_bytes, 4, n_threads, EM_LINES);
    return outfile;
}

string EM_delete_duplicate_lines_from_sorted_file(string infile){
    string outfile = get_temp_file_manager().create_filename();

    throwing_ifstream in(infile);
    throwing_ofstream out(outfile);

    LL line_count = 0;
    string line;
    string prev_line;
    while(in.getline(line)){
        if(line_count == 0 || line != prev_line)
            out << line << "\n";
        prev_line = line;
        line_count++;
    }

    return outfile;
}

string debug_binary_to_string_form(string infile){
    string outfile = get_temp_file_manager().create_filename();

    throwing_ifstream in(infile,ios::binary);
    throwing_ofstream out(outfile);

    char cur[8+8+8]; // size parameter and two long longs
    
    while(in.read(cur, 8+8+8)){
        //LL size = parse_big_endian_LL(cur);
        LL x = parse_big_endian_LL(cur+8);
        LL y = parse_big_endian_LL(cur+16);
        out << x << ";" << y << "\n";
    }

    return outfile;
}

string EM_delete_duplicate_LL_pair_records(string infile){
    string outfile = get_temp_file_manager().create_filename();

    throwing_ifstream in(infile, ios::binary);
    throwing_ofstream out(outfile, ios::binary);

    char prev[8+8]; // two long longs
    char cur[8+8]; // two long longs
    
    LL record_count = 0;
    while(in.read(cur, 8+8)){

        if(record_count == 0 || memcmp(prev, cur, 8+8) != 0){
            // The first record or different from the previous record
            out.write(cur, 8+8);
        }

        memcpy(prev, cur, 8+8);
        record_count++;
    }

    return outfile;
}

// (node, color) pairs to (node, color list) pairs
string EM_collect_colorsets_binary(string infile){
    string outfile = get_temp_file_manager().create_filename();

    throwing_ifstream in(infile, ios::binary);
    throwing_ofstream out(outfile, ios::binary);

    LL active_key = -1;
    vector<LL> cur_value_list;

    char buffer[8+8];

    while(true){
        in.read(buffer,8+8);
        if(in.stream.eof()) break;

        LL key = parse_big_endian_LL(buffer + 0);
        LL value = parse_big_endian_LL(buffer + 8);

        if(key == active_key) cur_value_list.push_back(value);
        else{
            if(active_key != -1){
                sort(cur_value_list.begin(), cur_value_list.end());
                LL record_size = 8 * (1 + 1 + cur_value_list.size());
                write_big_endian_LL(out, record_size);
                write_big_endian_LL(out, active_key);
                for(LL x : cur_value_list){
                    write_big_endian_LL(out, x);
                }
            }
            active_key = key;
            cur_value_list.clear();
            cur_value_list.push_back(value);
        }
    }

    // Last one
    if(active_key != -1){
        sort(cur_value_list.begin(), cur_value_list.end());
        LL record_size = 8 * (1 + 1 + cur_value_list.size());
        write_big_endian_LL(out, record_size);
        write_big_endian_LL(out, active_key);
        for(LL x : cur_value_list){
            write_big_endian_LL(out, x);
        }
    }

    return outfile;
}

// Input: records (node, color set), in the format where
// first we have the length of the whole record in bytes
string EM_sort_by_colorsets_binary(string infile, LL ram_bytes, LL n_threads){

    auto cmp = [&](const char* x, const char* y) -> bool{
        LL nx = parse_big_endian_LL(x);
        LL ny = parse_big_endian_LL(y);
        LL c = memcmp(x + 8 + 8, y + 8 + 8, min(nx-8-8,ny-8-8));
        if(c < 0) return true;
        else if(c > 0) return false;
        else { // c == 0
            return nx < ny;
        }
    };
    
    string outfile = get_temp_file_manager().create_filename();
    EM_sort(infile, outfile, cmp, ram_bytes, 4, n_threads, EM_VARIABLE_BINARY);
    return outfile;
}

string EM_collect_nodes_by_colorset_binary(string infile){
    string outfile = get_temp_file_manager().create_filename();

    throwing_ifstream in(infile, ios::binary);
    throwing_ofstream out(outfile, ios::binary);

    vector<LL> active_key;
    vector<LL> cur_value_list;

    vector<char> buffer(8);
    vector<LL> key; // Reusable space
    while(true){
        key.clear();
        in.read(buffer.data(),8);
        if(in.stream.eof()) break;
        LL record_len = parse_big_endian_LL(buffer.data());
        assert(record_len >= 8+8+8);
        while(buffer.size() < record_len)
            buffer.resize(buffer.size()*2);
        in.read(buffer.data()+8,record_len-8); // Read the rest
        LL value = parse_big_endian_LL(buffer.data()+8);
        LL color_set_size = (record_len - 8 - 8) / 8;
        
        for(LL i = 0; i < color_set_size; i++){
            key.push_back(parse_big_endian_LL(buffer.data() + 8 + 8 + i*8));
        }
        if(key == active_key) cur_value_list.push_back(value);
        else{
            if(active_key.size() != 0){
                sort(cur_value_list.begin(), cur_value_list.end());
                LL record_size = 8 + 8 + 8*cur_value_list.size() + 8*active_key.size();
                assert(record_size > 0);

                // record = (record length, number of nodes, node list, color list)
                write_big_endian_LL(out, record_size);
                write_big_endian_LL(out, cur_value_list.size());
                for(LL x : cur_value_list){
                    write_big_endian_LL(out, x);
                }
                for(LL x : active_key){
                    write_big_endian_LL(out, x);
                }
            }
            active_key = key;
            cur_value_list.clear();
            cur_value_list.push_back(value);
        }
    }

    // Last one
    if(active_key.size() != 0){
        sort(cur_value_list.begin(), cur_value_list.end());
        LL record_size = 8 + 8 + 8*cur_value_list.size() + 8*active_key.size();
        write_big_endian_LL(out, record_size);
        write_big_endian_LL(out, cur_value_list.size());
        for(LL x : cur_value_list){
            write_big_endian_LL(out, x);
        }
        for(LL x : active_key){
            write_big_endian_LL(out, x);
        }
    }

    return outfile;
}


// Input: a file with ';'-separated distinct string pairs (x,y),
//        sorted by x if key = 0, and by y if key = 1
// Output: 
// if key_pos == 0:
//   then pairs (x,Y), where Y is a space-separated list of y which appear with x
// if key_pos == 1:
//   then pairs (X,y), where X is a space-separated list of x which appear with y
string EM_collect(string infile, LL key_pos){
    assert(key_pos == 0 || key_pos == 1);

    string outfile = get_temp_file_manager().create_filename();

    throwing_ifstream in(infile);
    throwing_ofstream out(outfile);

    string line;
    string cur_key = "";
    vector<string> cur_value_list;

    while(in.getline(line)){
        vector<string> tokens = split(line,';');
        assert(tokens.size() == 2);
        string key, value;
        if(key_pos == 0){
            key = tokens[0]; value = tokens[1];
        } else{
            key = tokens[1]; value = tokens[0];
        }
        assert(key != ""); assert(value != "");
        if(key == cur_key) cur_value_list.push_back(value);
        else{
            if(cur_key != ""){
                if(key_pos == 0) out << cur_key << ";";
                for(LL i = 0; i < cur_value_list.size(); i++){
                    out << cur_value_list[i];
                    if(i != cur_value_list.size() - 1) out << " ";
                }
                if(key_pos == 1) out << ";" << cur_key;
                out << "\n";
            }
            cur_key = key;
            cur_value_list.clear();
            cur_value_list.push_back(value);
        }
    }

    // Last one
    if(cur_key != ""){
        if(key_pos == 0) out << cur_key << ";";
        for(LL i = 0; i < cur_value_list.size(); i++){
            out << cur_value_list[i];
            if(i != cur_value_list.size() - 1) out << " ";
        }
        if(key_pos == 1) out << ";" << cur_key;
        out << "\n";
    }

    return outfile;
}
