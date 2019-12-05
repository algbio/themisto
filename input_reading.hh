#pragma once
#include <utility>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "throwing_streams.hh"

using namespace std; // Bad practice but whatever

const int64_t FASTA_MODE = 0;
const int64_t FASTQ_MODE = 1;

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
    
    throwing_ifstream* file;

public:

    string header;
    int64_t mode;
    string dummy; // Used to skip over lines in fastq
    
    // mode is FASTA_MODE of FASTQ_MODE defined in this file
    Read_stream(throwing_ifstream* file, string header, int64_t mode) : file(file), header(header), mode(mode) {
    
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
            return true;
        } else if(mode == FASTQ_MODE){
            int next_char = file->stream.peek();
            if(next_char == '\n' || next_char == '\r') {
                // End of read. Rewind two lines forward to get to the header of the next read
                // for the next read stream
                file->getline(dummy); // Consume the newline
                file->getline(dummy); // Consume the '+'-line
                file->getline(dummy); // Consume the quality values
                return false;
            }
            else{
                file->read(&c,1);
                return true;
            }
        } else{
            cerr << "Invalid sequence read mode: " << mode << endl;
            assert(false);
        }
    }

    string get_all(){ // todo: make more efficient?
        char c;
        string read;
        while(getchar(c)) read += c;
        return read;
    }

};

class Sequence_Reader{

public:

    throwing_ifstream file;
    int64_t mode;

    // mode is FASTA_MODE of FASTQ_MODE defined in this file
    Sequence_Reader(string filename, int64_t mode) : file(filename, ios::in | ios::binary), mode(mode) {
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
        Read_stream rs(&file, header,mode);
        return rs;
    }

    bool done(){
        return file.stream.peek() == EOF;
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

