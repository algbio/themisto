#pragma once
#include <utility>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "throwing_streams.hh"

using namespace std; // Bad practice but whatever

class Raw_file_stream{
private:
    
    throwing_ifstream file;
    
public:
    
    Raw_file_stream(string filename) : file(filename, ios::in | ios::binary) {}
    
    bool getchar(char& c){
        if(file.stream.eof()) return false; // End of stream
        file.read(&c,1); // Read 1 byte
        if(file.stream.eof()) return false;
        if(c == '\n' || c == '\r')
            std::cerr << "Warning: file contains a newline character" << std::endl;
        return true;
    }
    
};

class Read_stream{
    
private:
    
    bool end_of_read;
    throwing_ifstream* file;

public:

    string header;
    
    Read_stream(throwing_ifstream* file, string header) : end_of_read(false), file(file), header(header) {
    
    }

    // Behaviour: Tries to read a byte to c. If eof or '>' was read,
    // return false and return false from here on. If c
    // is '\n' or '\r', read another byte recursively.
    bool getchar(char& c){
        start:
        int next_char = file->stream.peek();
        if(next_char == EOF || next_char == '>') return false;
        if(next_char == '\n' || next_char == '\r'){
            file->read(&c,1);
            goto start; // "recursive call"
        }

        file->read(&c,1);
        return true;
    }

    string get_all(){ // todo: make more efficient?
        char c;
        string read;
        while(getchar(c)) read += c;
        return read;
    }

};



class FASTA_reader{
    
private:
    
    throwing_ifstream file;
    
public:
    
    FASTA_reader(string filename) : file(filename, ios::in | ios::binary) {}
    
    bool done(){
        return file.stream.peek() == EOF;
    }
    
    Read_stream get_next_query_stream(){
        string header;
        file.getline(header);
        Read_stream rs(&file, header);
        return rs;
    }
};

// Vector of (read, header) pairs
std::vector<std::pair<std::string, std::string> > parse_FASTA(std::string filename){
    throwing_ifstream input(filename);

    std::vector<std::pair<std::string,std::string> > reads;

    std::string line;
    while(input.getline(line)){
        while(line.size() > 0 && isspace(line.back()))
            line.pop_back(); // Trim trailing whitespace just in case

        if(line.size() == 0 && !(input.stream.eof())) continue; // Ignore empty lines, just in case

        if(line[0] == '>')
            reads.push_back({"", line}); // Start new read
        else 
            reads.back().first += line; // Append base pairs to read
    }

    //Erase invalid reads
    //reads.erase(std::remove_if(reads.begin(), reads.end(), is_invalid_read), reads.end());
    return reads;
}

