#pragma once

/*
  Buffered reading for FASTA and FASTQ files.
  Authors: Jarno Alanko & Simon Puglisi
*/

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

struct ReadBatch{
    public:
    uint64_t firstReadID; //read id of first read in the batch
    vector<char> data; //concatenation of 
    vector<uint32_t> readStarts; //indicates starting positions of reads in data array, last item is a dummy
};

class ReadBatchReader{
public:
    virtual ReadBatch *getNextReadBatch() = 0;
    virtual ~ReadBatchReader() {}
};

class BufferedFastqStreamReader : public ReadBatchReader{

private:

    ifstream* file;
    unsigned char *buffer;
    uint64_t bufferCap; 
    uint64_t bufferSize;
    vector<char> readPrefix;
    uint64_t lineNumber; 
    uint64_t state;
    char prevchar;
    #define STATE_READING_HEADER 0
    #define STATE_READING_READ 1
    #define STATE_READING_PLUS 2
    #define STATE_READING_QS 3

public:

    // mode is FASTA_MODE of FASTQ_MODE defined in this file
    BufferedFastqStreamReader(ifstream* file, uint64_t bufferCap) : file(file), bufferCap(bufferCap){
       buffer = new unsigned char[bufferCap];
       lineNumber = 0;
       prevchar = 0;
    }

    ~BufferedFastqStreamReader(){
       delete [] buffer;
    }

    ReadBatch *getNextReadBatch(){
       if(file->eof()) return nullptr;
       file->read((char *)buffer,bufferCap);
       bufferSize = file->gcount(); 
       ReadBatch *rb = new ReadBatch();
       //cerr << "1.\n";
       //if an open read exists, copy its prefix onto the read batch
       if(readPrefix.size()){
          rb->readStarts.push_back(0);
          for(uint64_t i=0; i<readPrefix.size(); i++){
             rb->data.push_back(readPrefix[i]);
          }
          readPrefix.resize(0);
       }
       //cerr << "2.\n";
       //now process what's in the buffer
       for(uint64_t i=0; i<bufferSize; i++){
          if(buffer[i] == '\n' || buffer[i] == '\r'){
             lineNumber++;
          }
          else if((lineNumber%4) == STATE_READING_READ){
             if(prevchar == '\n' || prevchar == '\r'){
                uint64_t start = rb->data.size();
                rb->readStarts.push_back(start);
             }
             rb->data.push_back(buffer[i]); //base
          }
          prevchar = buffer[i];
       }
       //cerr << "3.\n";
       //all chars in buffer have been processed
       //but were we halfway through reading a read?
       if((lineNumber%4) == STATE_READING_READ){
          //cerr << "Open read, readPrefix.size(): " << readPrefix.size() << '\n';
          //do something appropriate
          for(uint64_t i = (rb->readStarts[rb->readStarts.size()-1]); i < (rb->data.size()); i++){
             readPrefix.push_back(rb->data[i]);
          }
       }else{
          //push a dummy onto readStarts
          rb->readStarts.push_back(rb->data.size());
       }
       //cerr << "4.\n";
       return rb;
    }

};

class ReadBatchIterator{
// An iterator reading one read starting from the given start position

public:

    ReadBatch* batch;
    uint64_t pos; // Current position

    ReadBatchIterator(ReadBatch* batch, int64_t start) : batch(batch), pos(start) {}

    // Returns (pointer to start, length). If done, return (NULL, 0)
    pair<const char*, uint64_t> getNextRead(){
        if(pos >= (batch->readStarts.size()-1)) return {NULL,0}; // End of batch
        uint64_t len = batch->readStarts[pos+1] - batch->readStarts[pos];
        uint64_t start = batch->readStarts[pos];
        pos++; //move to next read ready for next call to this function
        return {(batch->data.data()) + start, len};
    }
};
