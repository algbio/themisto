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
#include <memory>
#include "buffered_streams.hh"

namespace seq_io{

const std::vector<std::string> fasta_suffixes = {".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa"};
const std::vector<std::string> fastq_suffixes = {".fastq", ".fq"};

// Table mapping ascii values of characters to their reverse complements,
// lower-case to lower case, upper-case to upper-case. Non-ACGT characters
// are mapped to themselves.
static constexpr unsigned char rc_table[256] =
{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
56, 57, 58, 59, 60, 61, 62, 63, 64, 84, 66, 71, 68, 69, 70, 67, 72, 73,
74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 65, 85, 86, 87, 88, 89, 90, 91,
92, 93, 94, 95, 96, 116, 98, 103, 100, 101, 102, 99, 104, 105, 106, 107,
108, 109, 110, 111, 112, 113, 114, 115, 97, 117, 118, 119, 120, 121, 122,
123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137,
138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152,
153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182,
183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212,
213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227,
228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242,
243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255};

enum Format {FASTA, FASTQ};

struct FileFormat{
    Format format;
    bool gzipped;
    std::string extension; // Includes the possible .gz extension
};

inline void reverse_complement_c_string(char* S, int64_t len){
    std::reverse(S, S + len);
    for(int64_t i = 0; i < len; i++)
        S[i] = rc_table[(int)S[i]];
}

// Inlined to keep this a header-only library
inline FileFormat figure_out_file_format(std::string filename){
    Format fasta_or_fastq;
    bool gzipped = false;
    std::string extension;

    std::string gzip_suffix = "";
    if(filename.size() >= 3 && filename.substr(filename.size()-3) == ".gz"){
        filename = filename.substr(0, filename.size()-3); // Drop .gz
        gzipped = true;
    }

    for(int64_t i = (int64_t)filename.size()-1; i >= 0; i--){
        if(filename[i] == '.'){
            std::string ending = filename.substr(i);
            
            if(std::find(fasta_suffixes.begin(), fasta_suffixes.end(), ending) != fasta_suffixes.end()){
                fasta_or_fastq = FASTA;
                extension = ending + (gzipped ? ".gz" : "");
                return {fasta_or_fastq, gzipped, extension};
            }
            
            if(std::find(fastq_suffixes.begin(), fastq_suffixes.end(), ending) != fastq_suffixes.end()){
                fasta_or_fastq = FASTQ;
                extension = ending + (gzipped ? ".gz" : "");
                return {fasta_or_fastq, gzipped, extension};
            }

            throw(std::runtime_error("Unknown file format: " + filename + (gzipped ? ".gz" : "")));
        }
    }
    throw(std::runtime_error("Unknown file format: " + filename + (gzipped ? ".gz" : "")));
}

class NullStream : public std::ostream {
public:
  NullStream() : std::ostream(nullptr) {}
};

template <class T>
const NullStream &operator<<(NullStream &&os, const T &value) { 
  return os;
}


// Return the filename of the reverse-complemented file
template<typename reader_t, typename writer_t>
void create_reverse_complement_file(const std::string& infile, const std::string& outfile){
    seq_io::FileFormat fileformat = seq_io::figure_out_file_format(infile);

    reader_t sr(infile);
    writer_t sw(outfile);

    while(true) {
        int64_t len = sr.get_next_read_to_buffer();
        if(len == 0) break;

        // Reverse complement
        char* buf = sr.read_buf;
        std::reverse(buf, buf + len);
        for(int64_t i = 0; i < len; i++) buf[i] = rc_table[(int)buf[i]];

        sw.write_sequence(buf, len);
    }
}

// Creates a reverse-complement version of each file and return the filenames of the new files
template<typename reader_t, typename writer_t>
void create_reverse_complement_files(const std::vector<std::string>& infiles, const std::vector<std::string>& outfiles){
    assert(infiles.size() == outfiles.size());
    for(int64_t i = 0; i < infiles.size(); i++){
        create_reverse_complement_file<reader_t, writer_t>(infiles[i], outfiles[i]);
    }
}

template<typename ifstream_t = Buffered_ifstream<std::ifstream>> // The underlying file stream.
class Reader {

// The class is used like this:
// Sequence_Reader_Buffered sr;
// while(true) { 
//   int64_t len = sr.get_next_read_to_buffer();
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

std::unique_ptr<ifstream_t> stream;
int64_t mode;
int64_t read_buf_cap;
int64_t header_buf_cap;

bool reverse_complements = false; // Whether reverse complements are enabled
bool return_rc_next = false; // If reverse complements are enabled, this flag is used internally to manage the process
std::string filename;

std::vector<char> rc_buf; // Internal buffer for reverse complements

std::string new_read_buf, new_header_buf, new_plus_buf, new_quality_buf;
std::string fasta_read_concat_buf;

public:

    // These buffers are intended to be read from outside the class
    char* read_buf; // Stores a sequence read
    char* header_buf; // Stores the header of a read (without the '>' or '@')

    void read_first_char_and_sanity_check(){
        
        char c = 0; stream->get(&c);
        if(mode == FASTA && c != '>')
            throw std::runtime_error("ERROR: FASTA file " + filename + " does not start with '>'");
        if(mode == FASTQ && c != '@')
            throw std::runtime_error("ERROR: FASTQ file " + filename + " does not start with '@'");

        // This leaves the input stream pointer after the first character, but
        // get_next_read_to_buffer is written such that it's ok.
    }

    // mode should be FASTA_MODE or FASTQ_MODE
    // Note: FASTQ mode does not support multi-line FASTQ
    Reader(std::string filename, int64_t mode) : mode(mode), filename(filename) {
        stream = std::make_unique<ifstream_t>(filename, std::ios::binary);
        if(mode != FASTA && mode != FASTQ)
            throw std::invalid_argument("Unkown sequence format");
        
        read_buf_cap = 256;
        read_buf = (char*)malloc(read_buf_cap);

        header_buf_cap = 256;
        header_buf = (char*)malloc(header_buf_cap);

        read_first_char_and_sanity_check();
    }

    Reader(std::string filename) : filename(filename) {
        stream = std::make_unique<ifstream_t>(filename, std::ios::binary);
        seq_io::FileFormat fileformat = figure_out_file_format(filename);
        if(fileformat.format == FASTA) mode = FASTA;
        else if(fileformat.format == FASTQ) mode = FASTQ;
        else throw(std::runtime_error("Unknown file format: " + filename));

        read_buf_cap = 256;
        read_buf = (char*)malloc(read_buf_cap);

        header_buf_cap = 256;
        header_buf = (char*)malloc(header_buf_cap);

        read_first_char_and_sanity_check();
    }

    void enable_reverse_complements(){
        reverse_complements = true;
        return_rc_next = false;
    }


    ~Reader(){
        free(read_buf);
        free(header_buf);
    }

    void rewind_to_start(){
        // Create a new stream
        stream = std::make_unique<ifstream_t>(filename, std::ios::binary);

        read_first_char_and_sanity_check();
        return_rc_next = false;
    }

    void grow_buf_if_needed(char** buf, int64_t* cap, int64_t required_size){
        while(*cap < required_size){
            *cap *= 2;
            *buf = (char*)realloc(*buf, *cap);
        }
    }

    int64_t get_mode() const {return mode;}

    // Returns length of read, or zero if no more reads.
    // The read is null-terminated.
    // The read is stored in the member pointer `read_buffer`
    // The header is stored in the member pointer `header buffer`
    // When called, the read that is currently in the buffer is overwritten
    int64_t get_next_read_to_buffer() {

        if(reverse_complements){
            if(return_rc_next){
                strcpy(read_buf, rc_buf.data());
                rc_buf.clear();
                return_rc_next = false;
                return strlen(read_buf);
            } else {
                if(!stream->eof()) return_rc_next = true;
            }
        }

        if(stream->eof()) return 0;

        if(mode == FASTA){
            fasta_read_concat_buf.clear();
            char c = 0;
            stream->getline(new_header_buf);
            if (stream->eof()) throw std::runtime_error("FASTA file " + filename + " ended unexpectedly.");

            // Read the sequence

            stream->get(&c);
            if(c == '\n') 
                throw std::runtime_error("Empty line in FASTA file " + filename + ".");
            else if(c == '>')
                throw std::runtime_error("Empty sequence in FASTA file " + filename + ".");

            while(c != '>'){
                fasta_read_concat_buf.push_back(c);
                stream->getline(new_read_buf);
                if (stream->eof()) throw std::runtime_error("FASTA file " + filename + " ended unexpectedly.");
                fasta_read_concat_buf.append(new_read_buf);
                stream->get(&c); // Start of the next line
                if(stream->eof()) break;
            }

            for(char& c : fasta_read_concat_buf) c = toupper(c);

            int64_t read_len = fasta_read_concat_buf.size();

            grow_buf_if_needed(&header_buf, &header_buf_cap, new_header_buf.size() + 1); // +1: null terminator
            memcpy(header_buf, new_header_buf.data(), new_header_buf.size());
            header_buf[new_header_buf.size()] = '\0';

            grow_buf_if_needed(&read_buf, &read_buf_cap, fasta_read_concat_buf.size() + 1); // +1: null terminator.
            memcpy(read_buf, fasta_read_concat_buf.data(), fasta_read_concat_buf.size());
            read_buf[fasta_read_concat_buf.size()] = '\0';

            if(reverse_complements){
                // Store the reverse complement for later
                for(int64_t i = 0; i < read_len+1; i++) // +1: also copy the null
                    rc_buf.push_back(read_buf[i]);
                reverse_complement_c_string(rc_buf.data(), read_len);
            }

            return read_len;
        } else if(mode == FASTQ){
            stream->getline(new_header_buf);
            if (stream->eof()) throw std::runtime_error("FASTQ file " + filename + " ended unexpectedly.");

            stream->getline(new_read_buf);
            if (stream->eof()) throw std::runtime_error("FASTQ file " + filename + " ended unexpectedly.");

            stream->getline(new_plus_buf);
            if (stream->eof()) throw std::runtime_error("FASTQ file " + filename + " ended unexpectedly.");

            stream->getline(new_quality_buf);
            if (stream->eof()) throw std::runtime_error("FASTQ file " + filename + " ended unexpectedly.");

            for(char& c : new_read_buf) c = toupper(c);

            grow_buf_if_needed(&header_buf, &header_buf_cap, new_header_buf.size() + 1); // +1: null terminator
            memcpy(header_buf, new_header_buf.data(), new_header_buf.size());
            header_buf[new_header_buf.size()] = '\0';

            grow_buf_if_needed(&read_buf, &read_buf_cap, new_read_buf.size() + 1); // +1: null terminator
            memcpy(read_buf, new_read_buf.data(), new_read_buf.size());
            read_buf[new_read_buf.size()] = '\0';

            int64_t read_len = new_read_buf.size();

            char c;
            stream->get(&c); // Consume the '@' of the next read. If no more reads left, sets the eof flag.
            if(read_len == 0) throw std::runtime_error("Error: empty sequence in FASTQ file.");

            if(reverse_complements){
                // Store the reverse complement for later
                for(int64_t i = 0; i < read_len+1; i++) // +1: also copy the null
                    rc_buf.push_back(read_buf[i]);
                reverse_complement_c_string(rc_buf.data(), read_len);
            }

            return read_len;
        } else{
            throw std::runtime_error("Should not come to this else-branch");
        }
    }

    // Slow
    std::string get_next_read(){
        int64_t len = get_next_read_to_buffer();
        std::string read = (len > 0 ? std::string(read_buf) : "");
        return read;
    }

};

// Produces reads from multiple files like it was a single file
template<typename reader_t = seq_io::Reader<>>
class Multi_File_Reader{

    public:

    char* read_buf; // Does not own this memory
    char* header_buf; // Does not own this memory

    std::vector<std::string> filenames;
    int64_t current_file_idx;
    std::unique_ptr<reader_t> reader;
    bool reverse_complements = false;

    Multi_File_Reader(const std::vector<std::string>& filenames) : filenames(filenames), current_file_idx(0){
        if(filenames.size() > 0){
            reader = std::make_unique<reader_t>(filenames[0]);
        }
    }

    int64_t get_next_read_to_buffer(){
        if(current_file_idx == filenames.size()) return 0; // All files processed

        int64_t len = reader->get_next_read_to_buffer();
        while(len == 0){ // End of file
            current_file_idx++;
            if(current_file_idx == filenames.size()) return 0; // All files processed
            reader = std::make_unique<reader_t>(filenames[current_file_idx]);
            if(reverse_complements) reader->enable_reverse_complements();
            len = reader->get_next_read_to_buffer();
        }

        this->read_buf = reader->read_buf; // Update pointer in case there was a realloc
        this->header_buf = reader->header_buf; // Update pointer in case there was a realloc
        return len;
    }

    void enable_reverse_complements(){
        reverse_complements = true;
        if(filenames.size() > 0) reader->enable_reverse_complements();
    }

    void rewind_to_start(){
        current_file_idx = 0;
        if(filenames.size() > 0){
            reader = std::make_unique<reader_t>(filenames[0]);
            if(reverse_complements) reader->enable_reverse_complements();
        }
    }
};


template<typename ofstream_t = Buffered_ofstream<std::ofstream>> // The underlying file stream.
class Writer{

    std::string fasta_header = ">\n";
    std::string fastq_header = "@\n";
    std::string newline = "\n";
    std::string plus = "+";

    public:

    ofstream_t out;
    int64_t mode;

    // Tries to figure out the format based on the file extension.
    Writer(std::string filename) : out(filename) {
        seq_io::FileFormat fileformat = figure_out_file_format(filename);
        if(fileformat.format == FASTA) mode = FASTA;
        else if(fileformat.format == FASTQ) mode = FASTQ;
        else throw(std::runtime_error("Unknown file format: " + filename));
    }

    void write_sequence(const char* seq, int64_t len){
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

inline int64_t count_sequences(const std::string& filename){
    int64_t count = 0;
    if(figure_out_file_format(filename).gzipped){
        Reader<Buffered_ifstream<zstr::ifstream>> reader(filename);
        while(reader.get_next_read_to_buffer()) count++;
    } else{
        Reader<Buffered_ifstream<std::ifstream>> reader(filename);
        while(reader.get_next_read_to_buffer()) count++;
    }
    return count;
}



} // namespace seq_io
