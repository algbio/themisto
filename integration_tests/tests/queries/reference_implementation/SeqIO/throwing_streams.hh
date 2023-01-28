#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <string>
#include <sstream>
#include <vector>

using namespace std;

class throwing_ofstream{

public:

    // No copying, just moving
    throwing_ofstream(throwing_ofstream const& other) = delete;
    throwing_ofstream& operator=(throwing_ofstream const& other) = delete;
    throwing_ofstream(throwing_ofstream&&) = default;
    throwing_ofstream& operator = (throwing_ofstream&&) = default;

    throwing_ofstream() {}

    string filename;
    ofstream stream;

    // If file open fails, throws an error
    throwing_ofstream(string filename, ios_base::openmode mode = ios_base::out) : filename(filename) {
        stream.open(filename, mode);
        if(!stream.good()){
            throw runtime_error("Error opening file: " + filename);
        }
    }

    // Throws an error if write failed
    void write(const char* data, int64_t n){
        stream.write(data, n);
        if(!stream.good()){
            throw(runtime_error("Error writing to file " + filename));
        }
    }

    // Throws an error if file open failed
    void open(string filename, ios_base::openmode mode = ios_base::out){
        stream.open(filename, mode);
        if(!stream.good()){
            throw runtime_error("Error opening file: " + filename);
        }
    }

    void close(){
        stream.close();
    }

    void flush(){
        stream.flush();
    }
};

template <typename T>
throwing_ofstream& operator<<(throwing_ofstream& os, const T& t){
    os.stream << t;
    if(!os.stream.good()) 
        throw runtime_error("Error writing type " + string(typeid(T).name()) + " to file " + os.filename);
    return os;
}

class throwing_ifstream{

    public:

    // No copying, just moving
    throwing_ifstream(throwing_ifstream const& other) = delete;
    throwing_ifstream& operator=(throwing_ifstream const& other) = delete;
    throwing_ifstream(throwing_ifstream&&) = default;
    throwing_ifstream& operator = (throwing_ifstream&&) = default;

    string filename;
    ifstream stream;

    throwing_ifstream() {}
    // If file open failed, throws an error
    throwing_ifstream(string filename, ios_base::openmode mode = ios_base::in) : filename(filename) {
        stream.open(filename, mode);
        if(!stream.good()){
            throw runtime_error("Error opening file: '" + filename + "'");
        }
    }

    // If read succeeded, returns true
    // If read failed due to EOF, returns false.
    // If read failed due to other problems, throws an error.
    bool getline(string& line){
        std::getline(stream, line);
        if(!stream.good() && !stream.eof()) 
            throw runtime_error("Error reading from file '" + filename + "'");

        return !stream.eof();
    }

    // If read succeeded, returns true
    // If read failed due to EOF, returns false. (this means we can use the pattern while(in.read(data,n)){...}
    // If read failed due to other problems, throws an error.
    bool read(char* data, int64_t n){
        stream.read(data, n);
        if(!stream.good() && !stream.eof()) 
            throw runtime_error("Error reading from file '" + filename + "'");
        return !stream.eof();
    }

    // Returns number of characters read by the previous read operation
    int64_t gcount(){
        return stream.gcount();
    }

    // If open fails, throws an error
    void open(string filename, ios_base::openmode mode = ios_base::in){
        stream.open(filename, mode);
        if(!stream.good()){
            throw runtime_error("Error opening file: '" + filename + "'");
        }
    }

    void close(){
        stream.close();
    }

};

// Returns the stream object
// If read failed and it's not due to EOF, throws an error
template <typename T>
throwing_ifstream& operator>>(throwing_ifstream& is, T& t){
    is.stream >> t;
    if(!is.stream.good() && !is.stream.eof()) 
        throw runtime_error("Error reading type " + string(typeid(T).name()) + " from file " + is.filename);
    return is;
}
