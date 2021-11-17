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

    // Return true is get(c) has returned false
    bool eof(){
        return is_eof;
    }

};

class Sequence_Reader_Buffered {

// The class is used like this:
// Sequence_Reader_Buffered sr;
// while(true) { 
//   len = get_next_read_to_buffer(buf)
//   if(len == 0) break;
//   do something with sr.read_buffer
//}
//
// or (slow):
// while(true) { 
//   read = sr.get_next_read()
//   if(read.size() == 0) break;
//}

BufferedStream stream;
LL mode;
LL read_buf_cap ;

public:

    char* read_buffer;

    // mode should be FASTA_MODE or FASTQ_MODE
    // Note: FASTQ mode does not support multi-line FASTQ
    Sequence_Reader_Buffered(string filename, LL mode) : stream(filename), mode(mode) {
        // todo: check that fasta files start with > and fastq files start with @
        if(mode != FASTA_MODE && mode != FASTQ_MODE)
            throw std::invalid_argument("Unkown sequence format");
        
        read_buf_cap = 256;
        read_buffer = (char*)malloc(read_buf_cap);
    }

    ~Sequence_Reader_Buffered(){
        free(read_buffer);
    }

    // Returns length of read, or zero if no more reads.
    // The read is null-terminated.
    // The read is stored in the member pointer `read_buffer`
    // When called, the read that is currently in the buffer is overwritten
    LL get_next_read_to_buffer() {
        
        if(stream.eof()){
            read_buffer = nullptr;
            return 0;
        }

        if(mode == FASTA_MODE){
            char c = 0;
            while(c != '\n') stream.get(&c); // Skip fasta header line

            LL buf_pos = 0;
            while(true){
                if(!stream.get(&c)) break;
                else {
                    if(c == '\n') continue;
                    else if(c == '>') break;
                    else {
                        if(buf_pos + 1 >= read_buf_cap) { // +1: space for null terminator
                            read_buf_cap *= 2;
                            read_buffer = (char*)realloc(read_buffer, read_buf_cap);
                        }
                        read_buffer[buf_pos++] = toupper(c);
                    }
                }
            }
            read_buffer[buf_pos] = '\0';
            return buf_pos;
        } else if(mode == FASTQ_MODE){
            char c = stream.get(&c);
            if(stream.eof()){
                read_buffer = nullptr;
                return 0;
            }
            while(c != '\n') stream.get(&c); // Skip header line
            LL buf_pos = 0;
            while(true){
                stream.get(&c);
                if(c == '\n') break; // End of read
                if(buf_pos + 1 >= read_buf_cap) { // +1: space for null terminator
                    read_buf_cap *= 2;
                    read_buffer = (char*)realloc(read_buffer, read_buf_cap);
                }
                read_buffer[buf_pos++] = toupper(c);
            }
            read_buffer[buf_pos] = '\0';

            c = 0;
            while(c != '\n') stream.get(&c); // Skip '+'-line

            c = 0;
            while(c != '\n') stream.get(&c); // Skip quality line

            return buf_pos;
        } else{
            throw std::runtime_error("Should not come to this else-branch");
        }
    }

    // Slow
    string get_next_read(){
        LL len = get_next_read_to_buffer();
        string read = (len > 0 ? string(read_buffer) : "");
        return read;
    }

};