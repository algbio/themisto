#pragma once

#include <fstream>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>

#include "zstr/zstr.hpp"

namespace seq_io{

// The c++ ifstream and ofstream classes are buffered. But each read involves a virtual function
// call, which can be slow if the reads or writes are in small chunks. The buffer is also pretty
// small by default. These classes store the input/output in a large buffer and call the stream
// only when the buffer is full, which results in far less virtual function calls and better performance.
// These classes take as template parameter the underlying ifstream or ofstream, so you can use any
// stream that has the same interface as the std streams.

template<class ifstream_t = std::ifstream> // The underlying ifstream
class Buffered_ifstream{

private:

    Buffered_ifstream(const Buffered_ifstream& temp_obj) = delete; // No copying
    Buffered_ifstream& operator=(const Buffered_ifstream& temp_obj) = delete;  // No copying


    std::vector<char> buf;
    int64_t buf_cap = 1 << 20;

    int64_t buf_pos = 0;
    int64_t buf_size = 0;
    bool is_eof = false;
    ifstream_t* inner = nullptr;

public:

    Buffered_ifstream(Buffered_ifstream&&) = default; // Movable
    Buffered_ifstream& operator = (Buffered_ifstream&&) = default;  // Movable

    Buffered_ifstream() {}
    ~Buffered_ifstream() {
        delete inner;
    }

    Buffered_ifstream(std::string filename, std::ios_base::openmode mode = std::ios_base::in) {
        inner = new ifstream_t(filename, mode);
        if(!inner->good()) throw std::runtime_error("Error opening file " + filename);
        buf.resize(buf_cap);
    }

    // Reads one byte to the given location.
    // Returns true if read was succesful (no eof)
    bool get(char* c){
        if(is_eof) return false;

        if(buf_pos < buf_size){
            *c = buf[buf_pos++];
        } else{
            // Try to refill buffer
            if (inner->eof()) {
                is_eof = true;
                return false;
            }
            inner->read(buf.data(), buf_cap);
            buf_size = inner->gcount();
            if(buf_size == 0){
                is_eof = true;
                return false;
            }

            buf_pos = 0;
            *c = buf[buf_pos++];
            return true;
        }
        return true;
    }

    // Read up to n bytes to dest and returns the number of bytes read
    int64_t read(char* dest, int64_t n){
        char* ptr = dest;
        for(int64_t i = 0; i < n; i++){
            if(!get(ptr)) break;
            ptr++;
        }
        return ptr - dest;
    }

    // Should work exactly like std::ifstream::getline
    bool getline(std::string& line){
        line.clear();
        while(true){
            char c; get(&c);
            if(eof()) return line.size() > 0; // Last line can end without a newline
            if(c == '\n') return true; // Do not push the newline but just return
            line.push_back(c);
        }
    }

    // Return true if get(c) has returned false
    bool eof(){
        return is_eof;
    }

    void open(std::string filename, std::ios_base::openmode mode = std::ios_base::in){
        delete inner;
        inner = new ifstream_t(filename, mode);
        buf.resize(buf_cap);
        if(!inner->good()) throw std::runtime_error("Error opening file " + filename);
        buf_size = 0;
        buf_pos = 0;
        is_eof = false;
    }

    void close(){
        delete inner;
        inner = nullptr;
    }


    void set_buffer_capacity(int64_t cap){
        this->buf_cap = cap;
        buf.resize(buf_cap);
    }

    void rewind_to_start(){
        inner->clear();
        inner->seekg(0);
        buf_pos = 0;
        buf_size = 0;
        is_eof = false;
    }

    void seekg(int64_t pos){
        buf_pos = 0;
        buf_size = 0;
        is_eof = false;

        inner->clear();
        inner->seekg(pos);
    }


};

template<class ofstream_t = std::ofstream> // The underlying ifstream
class Buffered_ofstream{

private:

    Buffered_ofstream(const Buffered_ofstream& temp_obj) = delete; // No copying
    Buffered_ofstream& operator=(const Buffered_ofstream& temp_obj) = delete;  // No copying

    std::vector<char> buf;
    int64_t buf_size = 0;
    int64_t buf_cap = 1 << 20;
    ofstream_t* stream = nullptr;

    void empty_internal_buffer_to_stream(){
        if(stream){
            stream->write(buf.data(), buf_size);
            buf_size = 0;
            if (!stream->good()) {
                throw std::runtime_error("Error writing to file");
            }
        }
    }

public:

    Buffered_ofstream(Buffered_ofstream&&) = default; // Movable
    Buffered_ofstream& operator = (Buffered_ofstream&&) = default;  // Movable

    Buffered_ofstream(){}

    Buffered_ofstream(std::string filename, std::ios_base::openmode mode = std::ios_base::out){
        stream = new ofstream_t(filename, mode);
        if(!stream->good()) throw std::runtime_error("Error opening file " + filename);
        buf.resize(buf_cap);
    }

    void set_buffer_capacity(int64_t cap){
        this->buf_cap = cap;
        buf.resize(buf_cap);
    }

    void write(const char* data, int64_t n){
        for(int64_t i = 0; i < n; i++){
            if(buf_cap == buf_size) empty_internal_buffer_to_stream();
            buf[buf_size++] = data[i];
        }
    }

    void open(std::string filename, std::ios_base::openmode mode = std::ios_base::out){
        delete stream;
        stream = new ofstream_t(filename, mode);
        if(!stream->good()) throw std::runtime_error("Error opening file " + filename);
        buf.resize(buf_cap);
        buf_size = 0;
    }

    void close(){
        empty_internal_buffer_to_stream();
        delete stream; // Flushes also
        stream = nullptr;
    }

    // Flush the internal buffer AND the file stream
    void flush(){
        empty_internal_buffer_to_stream();
        stream->flush();
    }

    ~Buffered_ofstream(){
        empty_internal_buffer_to_stream();
        delete stream;
    }

};

}
