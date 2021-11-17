#pragma once

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include "globals.hh"

using namespace std;

class BufferedStream{

    static const LL buf_cap = (1 << 20);
    vector<char> buf;
    LL buf_pos = 0;
    LL buf_size = 0;
    bool is_eof = false;
    ifstream stream;

    public:

    BufferedStream(string filename) : stream(filename){
        if(!stream.good()) throw std::runtime_error("Error opening file " + filename);
        buf.resize(buf_cap);
    }

    // returns true if read was succesful
    // Stores the read character into the given pointer location
    bool get(char* c){
        if(is_eof) return false;
        if(buf_pos == buf_size){
            stream.read(buf.data(), buf_cap);
            buf_size = stream.gcount();
            buf_pos = 0;
            if(buf_size == 0){
                is_eof = true;
                return false;
            }
        }
        *c = buf[buf_pos++];
        return true;
    }

    bool eof(){
        return is_eof;
    }

};

class Sequence_Reader_New {

BufferedStream stream;

public:
    Sequence_Reader_New(string fastafile) : stream(fastafile) {}

    // Returns length of read, or zero if no more reads.
    // The read is null-terminated.
    LL get_next_read_to_buffer(char** buffer, LL* buf_cap) {
        // reallocs buffer if needed

        assert(*buf_cap > 0);

        if(stream.eof()){
            *buffer = nullptr;
            return 0;
        }

        char c = 0;
        while(c != '\n') stream.get(&c); // Skip fasta header line

        LL buf_pos = 0;
        while(true){
            if(!stream.get(&c)) break;
            else {
                if(c == '\n') continue;
                else if(c == '>') break;
                else {
                    if(buf_pos + 1 >= *buf_cap) { // +1: space for null terminator
                        *buf_cap = *buf_cap * 2;
                        *buffer = (char*)realloc(*buffer, *buf_cap);
                    }
                    (*buffer)[buf_pos++] = toupper(c);
                }
            }
        }
        (*buffer)[buf_pos] = '\0';
        return buf_pos;

    }

    // Slow
    string get_next_read(){
        LL buf_cap = 256;
        char* buffer = (char*)malloc(buf_cap);
        LL len = get_next_read_to_buffer(&buffer, &buf_cap);
        string read = (len > 0 ? string(buffer) : "");
        free(buffer);
        return read;
    }

    // There is no done(). The class is called like:
    // while(true) { 
    //   len = get_next_read_to_buffer(buf)
    //   if(len == 0) break;
    //}
    // or
    // while(true) { 
    //   read = get_next_read_to_buffer(buf)
    //   if(read.size() == 0) break;
    //}
};