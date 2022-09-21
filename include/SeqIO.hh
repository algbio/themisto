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
#include <algorithm>
#include "globals.hh"
#include "throwing_streams.hh"
#include "buffered_streams.hh"

using namespace std;

namespace SeqIO{

enum Format {FASTA, FASTQ};

struct FileFormat{
    Format format;
    bool gzipped;
    string extension; // Includes the possible .gz extension
};

FileFormat figure_out_file_format(string filename);

class NullStream : public std::ostream {
public:
  NullStream() : std::ostream(nullptr) {}
};

template <class T>
const NullStream &operator<<(NullStream &&os, const T &value) { 
  return os;
}

// Creates a reverse-complement version of each file and return the filenames of the new files
template<typename reader_t, typename writer_t>
vector<string> create_reverse_complement_files(const vector<string>& files){
    vector<string> newfiles;
    for(string f : files){
        SeqIO::FileFormat fileformat = SeqIO::figure_out_file_format(f);

        string f_rev = get_temp_file_manager().create_filename("", fileformat.extension);
        newfiles.push_back(f_rev);

        reader_t sr(f);
        writer_t sw(f_rev);

        while(true) {
            LL len = sr.get_next_read_to_buffer();
            if(len == 0) break;

            // Reverse complement
            char* buf = sr.read_buf;
            std::reverse(buf, buf + len);
            for(LL i = 0; i < len; i++) buf[i] = get_rc(buf[i]);

            sw.write_sequence(buf, len);
        }
    }
    return newfiles;
}

template<typename ifstream_t = Buffered_ifstream<std::ifstream>> // The underlying file stream.
class Reader {

// The class is used like this:
// Sequence_Reader_Buffered sr;
// while(true) { 
//   LL len = sr.get_next_read_to_buffer();
//   if(len == 0) break;
//   do something with sr.read_buf
//}
//
// or (slow):
// while(true) { 
//   read = sr.get_next_read()
//   if(read.size() == 0) break;
//}

private:

Reader(const Reader& temp_obj) = delete; // No copying
Reader& operator=(const Reader& temp_obj) = delete;  // No copying

ifstream_t stream;
LL mode;
LL read_buf_cap;

public:

    char* read_buf;

    void read_first_char_and_sanity_check(){
        
        char c = 0; stream.get(&c);
        if(mode == FASTA && c != '>')
            throw runtime_error("ERROR: FASTA file does not start with '>'");
        if(mode == FASTQ && c != '@')
            throw runtime_error("ERROR: FASTQ file does not start with '@'");

        // This leaves the input stream pointer after the first character, but
        // get_next_read_to_buffer is written such that it's ok.
    }

    // mode should be FASTA_MODE or FASTQ_MODE
    // Note: FASTQ mode does not support multi-line FASTQ
    Reader(string filename, LL mode) : stream(filename, ios::binary), mode(mode) {
        // todo: check that fasta files start with > and fastq files start with @
        if(mode != FASTA && mode != FASTQ)
            throw std::invalid_argument("Unkown sequence format");
        
        read_buf_cap = 256;
        read_buf = (char*)malloc(read_buf_cap);

        read_first_char_and_sanity_check();
    }

    // mode should be FASTA_MODE or FASTQ_MODE
    // Note: FASTQ mode does not support multi-line FASTQ
    Reader(string filename) : stream(filename, ios::binary) {
        SeqIO::FileFormat fileformat = figure_out_file_format(filename);
        if(fileformat.format == FASTA) mode = FASTA;
        else if(fileformat.format == FASTQ) mode = FASTQ;
        else throw(runtime_error("Unknown file format: " + filename));

        // todo: check that fasta files start with > and fastq files start with @
        read_buf_cap = 256;
        read_buf = (char*)malloc(read_buf_cap);

        read_first_char_and_sanity_check();
    }


    ~Reader(){
        free(read_buf);
    }

    LL get_mode() const {return mode;}

    // Returns length of read, or zero if no more reads.
    // The read is null-terminated.
    // The read is stored in the member pointer `read_buffer`
    // When called, the read that is currently in the buffer is overwritten
    LL get_next_read_to_buffer() {
        
        if(stream.eof()){
            read_buf = nullptr;
            return 0;
        }

        if(mode == FASTA){
            char c = 0;
            while(c != '\n') stream.get(&c); // Skip fasta header line

            LL buf_pos = 0;
            while(true){
                if(!stream.get(&c)) break; // Last read end
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
            if(buf_pos == 0) throw std::runtime_error("Error: empty sequence in FASTA file.");
            read_buf[buf_pos] = '\0';
            return buf_pos;
        } else if(mode == FASTQ){
            char c = 0; stream.get(&c);
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

            if(buf_pos == 0) throw std::runtime_error("Error: empty sequence in FASTQ file.");
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

template<typename ofstream_t = Buffered_ofstream<std::ofstream>> // The underlying file stream.
class Writer{

    string fasta_header = ">\n";
    string fastq_header = "@\n";
    string newline = "\n";
    string plus = "+";

    public:

    ofstream_t out;
    LL mode;

    // Tries to figure out the format based on the file extension.
    Writer(string filename) : out(filename) {
        SeqIO::FileFormat fileformat = figure_out_file_format(filename);
        if(fileformat.format == FASTA) mode = FASTA;
        else if(fileformat.format == FASTQ) mode = FASTQ;
        else throw(runtime_error("Unknown file format: " + filename));
    }

    void write_sequence(const char* seq, LL len){
        if(mode == FASTA){
            // FASTA format
            out.write(fasta_header.c_str(), 2);
            out.write(seq, len);
            out.write(newline.c_str(), 1);
        } else{
            // FASTQ
            out.write(fastq_header.c_str(), 2);
            out.write(seq, len);
            out.write(newline.c_str(), 1);
            out.write(plus.c_str(), 1);
            out.write(newline.c_str(), 1);
            out.write(seq, len); // Use the read again for the quality values
            out.write(newline.c_str(), 1);
        }
    }

    // Flush the stream. The stream is also automatically flushed when the object is destroyed.
    void flush(){
        out.flush();
    }
};

/*

LEGACY UNBUFFERED INPUT READING BELOW.

*/

class Unbuffered_Read_stream{
    
private:
    
    throwing_ifstream* file;

public:

    string header;
    int64_t mode;
    string dummy; // Used to skip over lines in fastq
    bool upper_case_enabled;
    
    // mode is FASTA_MODE of FASTQ_MODE defined in this file
    Unbuffered_Read_stream(throwing_ifstream* file, string header, int64_t mode, bool upper_case_enabled) : file(file), header(header), mode(mode), upper_case_enabled(upper_case_enabled) {
    
    }

    // Behaviour: Tries to read a char to c by peeking the file stream. Return false
    // if could not get a character because the sequence ended, otherwise return true.
    // If this returns false then the file stream will be put in such a state that the
    // next character is the first character of the header of the next read (or the EOF
    // character if it was the last read).
    bool getchar(char& c){
        if(mode == FASTA){
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
        } else if(mode == FASTQ){
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

// Unbuffered!! If you don't need headers, use SeqIO::Reader
class Unbuffered_Reader{

public:

    throwing_ifstream file;
    int64_t mode;
    bool upper_case_enabled;

    void sanity_check(){
        if(mode == FASTA) {
            if(file.stream.peek() != '>'){
                throw runtime_error("Error: FASTA-file does not start with '>'");
            }
        }
        if(mode == FASTQ) {
            if(file.stream.peek() != '@'){
                throw runtime_error("Error: FASTQ-file does not start with '@'");
            }
        }
    }

    // mode is FASTA_MODE of FASTQ_MODE defined in this file
    Unbuffered_Reader(string filename, int64_t mode) : file(filename, ios::in | ios::binary), mode(mode), upper_case_enabled(true) {
        sanity_check();
    }

    Unbuffered_Reader(string filename) : file(filename, ios::in | ios::binary), upper_case_enabled(true) {
        SeqIO::FileFormat fileformat = figure_out_file_format(filename);
        if(fileformat.format == FASTA) mode = FASTA;
        else if(fileformat.format == FASTQ) mode = FASTQ;
        else throw(runtime_error("Unknown file format: " + filename));
        sanity_check();
    }

    Unbuffered_Read_stream get_next_query_stream(){
        string header;
        file.getline(header);
        if(header.size() < 1) throw runtime_error("Error: FASTA or FASTQ parsing: header does not start with '>' or '@'");
        header = header.substr(1); // Drop the '>' in FASTA or '@' in FASTQ
        Unbuffered_Read_stream rs(&file, header, mode, upper_case_enabled);
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

} // Namespace SeqIO
