#pragma once

#include "stdlib.h"
#include <iostream>
#include <string>
#include <random>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <thread>
#include <mutex>
#include <cstring>
#include <set>
#include <filesystem>

using namespace std;

class Temp_File_Manager{

// Cleans up all allocated files in the end automatically

private:

    string get_random_string(int64_t length){
        string S;
        for(int64_t i = 0; i < length; i++){
            S += alphabet[dist(urandom)];
        }
        return S;
    }

    std::uniform_int_distribution<int64_t> dist;
    std::random_device urandom;
    string temp_dir;
    vector<char> alphabet;
    set<string> used_names;
    std::mutex mutex;

    void check_dir_exists(string path){
        struct stat info;    
        if( stat( path.c_str(), &info ) != 0 ){
            printf( "Error: cannot access %s\n", path.c_str() );
            exit(1);
        }
        else if( info.st_mode & S_IFDIR ){
            // All good
        }    
        else{
            printf( "Error: %s is not a directory\n", path.c_str() );
            exit(1);
        }
    }

public:


    Temp_File_Manager() : urandom("/dev/urandom") {
        for(char c = 'a'; c <= 'z'; c++) alphabet.push_back(c);
        for(char c = 'A'; c <= 'Z'; c++) alphabet.push_back(c);
        for(char c = '0'; c <= '9'; c++) alphabet.push_back(c);
        dist = std::uniform_int_distribution<int64_t>(0, alphabet.size()-1);
    }

    void set_dir(string temp_dir){ // Needs to be called before calling get_temp_file_name
        check_dir_exists(temp_dir);
        this->temp_dir = temp_dir;
    }

    string get_dir(){
        return this->temp_dir;
    }

    string create_filename(){
        return create_filename("","");
    }

    string create_filename(string prefix){
        return create_filename(prefix, "");
    }

    string create_filename(string prefix, string suffix){
        // Make sure only one thread runs in this function at once
        std::lock_guard<std::mutex> lg(mutex);
        if(temp_dir == ""){
            cerr << "Error: temp dir not set" << endl;
            exit(1);
        }
        while(true){
            string name = temp_dir + "/" + prefix + get_random_string(25) + suffix; // 62^25 >= 10^44 different possibilities
            auto desc = open(name.c_str(), O_CREAT | O_EXCL, S_IRWXU); // Fails if file exists, otherwise creates it
            if(desc != -1){
                used_names.insert(name);
                close(desc);
                // Success. Let's delete the empty file that was created
                std::filesystem::remove(name);
                return name;
            } else if(errno != EEXIST){
                cerr << std::strerror(errno) << " " << name << endl;
                close(desc);
                exit(1);
            }
            
        }
    } // The lock guard goes out of scope and is destructed

    // delete_file: delete a file before the end of the lifetime of the manager
    void delete_file(string filename){
        if(used_names.count(filename) == 0){
            cerr << "Error: tried to delete a temp file that thet temp file manager did not create: " << filename << endl;
            exit(1);
        }
        std::filesystem::remove(filename.c_str());
        used_names.erase(filename);
    }

    void delete_all_files(){
        for(string name : used_names) remove(name.c_str());
        used_names.clear();
    }

    ~Temp_File_Manager(){
        delete_all_files();
    }

};