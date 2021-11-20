#pragma once

/*
  Buffered reading for FASTA and FASTQ files.
  Authors: Jarno Alanko & Simon Puglisi
*/

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include "globals.hh"
#include "throwing_streams.hh"
#include "buffered_streams.hh"

using namespace std;

const int64_t FASTA_MODE = 0;
const int64_t FASTQ_MODE = 1;

class NullStream : public std::ostream {
public:
  NullStream() : std::ostream(nullptr) {}
};

template <class T>
const NullStream &operator<<(NullStream &&os, const T &value) { 
  return os;
}

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

private:

Sequence_Reader_Buffered(const Sequence_Reader_Buffered& temp_obj) = delete; // No copying
Sequence_Reader_Buffered& operator=(const Sequence_Reader_Buffered& temp_obj) = delete;  // No copying

Buffered_ifstream stream;
LL mode;
LL read_buf_cap;

public:

    char* read_buf;

    // mode should be FASTA_MODE or FASTQ_MODE
    // Note: FASTQ mode does not support multi-line FASTQ
    Sequence_Reader_Buffered(string filename, LL mode) : stream(filename), mode(mode) {
        // todo: check that fasta files start with > and fastq files start with @
        if(mode != FASTA_MODE && mode != FASTQ_MODE)
            throw std::invalid_argument("Unkown sequence format");
        
        read_buf_cap = 256;
        read_buf = (char*)malloc(read_buf_cap);
    }

    ~Sequence_Reader_Buffered(){
        free(read_buf);
    }

    // Returns length of read, or zero if no more reads.
    // The read is null-terminated.
    // The read is stored in the member pointer `read_buffer`
    // When called, the read that is currently in the buffer is overwritten
    LL get_next_read_to_buffer() {
        
        if(stream.eof()){
            read_buf = nullptr;
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
                            read_buf = (char*)realloc(read_buf, read_buf_cap);
                        }
                        read_buf[buf_pos++] = toupper(c);
                    }
                }
            }
            read_buf[buf_pos] = '\0';
            return buf_pos;
        } else if(mode == FASTQ_MODE){
            char c = stream.get(&c);
            if(stream.eof()){
                read_buf = nullptr;
                return 0;
            }
            while(c != '\n') stream.get(&c); // Skip header line
            LL buf_pos = 0;
            while(true){
                stream.get(&c);
                if(c == '\n') break; // End of read
                if(buf_pos + 1 >= read_buf_cap) { // +1: space for null terminator
                    read_buf_cap *= 2;
                    read_buf = (char*)realloc(read_buf, read_buf_cap);
                }
                read_buf[buf_pos++] = toupper(c);
            }
            read_buf[buf_pos] = '\0';

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
        string read = (len > 0 ? string(read_buf) : "");
        return read;
    }

};

/*

LEGACY UNBUFFERED INPUT READING BELOW.

*/

class Read_stream{
    
private:
    
    throwing_ifstream* file;

public:

    string header;
    int64_t mode;
    string dummy; // Used to skip over lines in fastq
    bool upper_case_enabled;
    
    // mode is FASTA_MODE of FASTQ_MODE defined in this file
    Read_stream(throwing_ifstream* file, string header, int64_t mode, bool upper_case_enabled) : file(file), header(header), mode(mode), upper_case_enabled(upper_case_enabled) {
    
    }

    // Behaviour: Tries to read a char to c by peeking the file stream. Return false
    // if could not get a character because the sequence ended, otherwise return true.
    // If this returns false then the file stream will be put in such a state that the
    // next character is the first character of the header of the next read (or the EOF
    // character if it was the last read).
    bool getchar(char& c){
        if(mode == FASTA_MODE){
            start:
            int next_char = file->stream.peek();
            if(next_char == EOF || next_char == '>') return false;
            if(next_char == '\n' || next_char == '\r'){
                file->read(&c,1);
                goto start; // "recursive call"
            }
            file->read(&c,1);
            if(upper_case_enabled) c = toupper(c);
            return true;
        } else if(mode == FASTQ_MODE){
            int next_char = file->stream.peek();
            if(next_char == '\n' || next_char == '\r') {
                // End of read. Rewind two lines forward to get to the header of the next read
                // for the next read stream
                file->getline(dummy); // Consume the newline
                assert(file->stream.peek() == '+');
                file->getline(dummy); // Consume the '+'-line
                file->getline(dummy); // Consume the quality values
                return false;
            }
            else{
                file->read(&c,1);
                if(upper_case_enabled) c = toupper(c);
                return true;
            }
        } else{
            throw(std::runtime_error("Invalid sequence read mode: " + mode));
        }
    }

    string get_all(){ // todo: make more efficient?
        char c;
        string read;
        while(getchar(c)) read += c;
        return read;
    }

};

// Unbuffered!! If you don't need headers, use Squence_Reader_Buffered
class Sequence_Reader{

public:

    throwing_ifstream file;
    int64_t mode;
    bool upper_case_enabled;

    // mode is FASTA_MODE of FASTQ_MODE defined in this file
    Sequence_Reader(string filename, int64_t mode) : file(filename, ios::in | ios::binary), mode(mode), upper_case_enabled(true) {
        if(mode == FASTA_MODE) {
            if(file.stream.peek() != '>'){
                throw runtime_error("Error: FASTA-file does not start with '>'");
            }
        }
        if(mode == FASTQ_MODE) {
            if(file.stream.peek() != '@'){
                throw runtime_error("Error: FASTQ-file does not start with '@'");
            }
        }
    }

    Read_stream get_next_query_stream(){
        string header;
        file.getline(header);
        if(header.size() < 1) throw runtime_error("Error: FASTA or FASTQ parsing: header does not start with '>' or '@'");
        header = header.substr(1); // Drop the '>' in FASTA or '@' in FASTQ
        Read_stream rs(&file, header, mode, upper_case_enabled);
        return rs;
    }

    // If flag is true, then query streams will upper case all sequences (off by default)
    void set_upper_case(bool flag){
        upper_case_enabled = flag;
    }

    bool done(){
        return file.stream.peek() == EOF;
    }

};