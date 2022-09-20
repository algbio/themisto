#pragma once

#include "globals.hh"

Temp_File_Manager& get_temp_file_manager();

// The c++ ifstream and ofstream classes are buffered. But each read involves a virtual function
// call, which can be slow if the reads or writes are in small chunks. The buffer is also pretty
// small by default. These classes store the input/output in a large buffer and call the stream 
// only when the buffer is full, which results in far less virtual function calls and better performance.

class Buffered_ifstream{

private:

    Buffered_ifstream(const Buffered_ifstream& temp_obj) = delete; // No copying
    Buffered_ifstream& operator=(const Buffered_ifstream& temp_obj) = delete;  // No copying

    
    vector<char> buf;
    LL buf_cap = 1 << 20;

    LL buf_pos = 0;
    LL buf_size = 0;
    bool is_eof = false;
    ifstream stream;

public:

    Buffered_ifstream(Buffered_ifstream&&) = default; // Movable
    Buffered_ifstream& operator = (Buffered_ifstream&&) = default;  // Movable

    Buffered_ifstream() {}

    Buffered_ifstream(string filename, ios_base::openmode mode = ios_base::in) : stream(filename, mode) {
        if(!stream.good()) throw std::runtime_error("Error opening file " + filename);
        buf.resize(buf_cap);
    }

    // Reads one byte to the given location.
    // Returns true if read was succesful
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

    // Read up to n bytes to dest and returns the number of bytes read
    LL read(char* dest, LL n){
        char* ptr = dest;
        for(LL i = 0; i < n; i++){
            if(!get(ptr)) break;
            ptr++;
        }
        return ptr - dest;
    }

    bool getline(string& line){
        line.clear();
        while(true){
            char c; get(&c);
            if(eof()) return line.size() > 0;
            if(c == '\n') return true;
            line.push_back(c);
        }
    }

    // Return true if get(c) has returned false
    bool eof(){
        return is_eof;
    }

    void open(string filename, ios_base::openmode mode = ios_base::in){
        buf.resize(buf_cap);
        stream.open(filename, mode);
        if(!stream.good()) throw std::runtime_error("Error opening file " + filename);
        buf_size = 0;
        buf_pos = 0;
        is_eof = false;
    }

    void close(){
        stream.close();
    }


    void set_buffer_capacity(LL cap){
        this->buf_cap = cap;
        buf.resize(buf_cap);
    }

};


class Buffered_ofstream{

private:

    Buffered_ofstream(const Buffered_ofstream& temp_obj) = delete; // No copying
    Buffered_ofstream& operator=(const Buffered_ofstream& temp_obj) = delete;  // No copying

    vector<char> buf;
    LL buf_size = 0;
    LL buf_cap = 1 << 23;
    ofstream stream;

    void empty_internal_buffer_to_stream(){
        try{
            stream.write(buf.data(), buf_size);
        } catch(const std::ios_base::failure& e){
            cerr << "Error writing to disk (perhaps the disk is full?)" << endl;
            cerr << e.what() << endl;
            cerr << e.code() << endl;
            get_temp_file_manager().delete_all_files();
            exit(1);
        }
        buf_size = 0;
    }

public:

    Buffered_ofstream(Buffered_ofstream&&) = default; // Movable
    Buffered_ofstream& operator = (Buffered_ofstream&&) = default;  // Movable

    Buffered_ofstream(){}
    Buffered_ofstream(string filename, ios_base::openmode mode = ios_base::out) : stream(filename, mode){
        if(!stream.good()) throw std::runtime_error("Error opening file " + filename);
        stream.exceptions(stream.failbit | stream.badbit);
        buf.resize(buf_cap);
    }

    void set_buffer_capacity(LL cap){
        this->buf_cap = cap;
        buf.resize(buf_cap);
    }

    void write(const char* data, int64_t n){
        for(LL i = 0; i < n; i++){
            if(buf_cap == buf_size) empty_internal_buffer_to_stream();
            buf[buf_size++] = data[i];
        }
    }

    void open(string filename, ios_base::openmode mode = ios_base::out){
        buf.resize(buf_cap);
        stream.open(filename, mode);
        if(!stream.good()) throw std::runtime_error("Error opening file " + filename);
        buf_size = 0;
    }

    void close(){
        empty_internal_buffer_to_stream();
        stream.close();
    }

    // Flush the internal buffer AND the file stream
    void flush(){
        empty_internal_buffer_to_stream();
        stream.flush();
    }

    ~Buffered_ofstream(){
        empty_internal_buffer_to_stream();
    }

};
