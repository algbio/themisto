#pragma once

#include <vector>
#include <string>
#include "SeqIO/SeqIO.hh"

using namespace std;

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


vector<string> readlines(string filename){
    vector<string> lines;
    string line;
    throwing_ifstream in(filename);
    while(getline(in.stream,line)){
        lines.push_back(line);
    }
    return lines;
}

vector<int64_t> read_colorfile(string filename){
    vector<string> lines = readlines(filename);
    vector<int64_t> colors;
    for(const string& line : lines){
        colors.push_back(stoll(line));
    }
    return colors;
}

int64_t count_sequences(string filename){
    SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in(filename);
    int64_t count = 0;
    while(in.get_next_read_to_buffer() > 0) count++;
    return count;
}

// Assumes gzipped data
vector<string> read_sequences(vector<string> filenames){
    vector<string> seqs;
    for(const string& f : filenames){
        SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> in(f);
        while(true){
            int64_t len = in.get_next_read_to_buffer();
            if(len == 0) break;
            seqs.push_back(string(in.read_buf));
        }
    }
    return seqs;
}

// Returns true iff S[start .. start + k - 1] has only characters A,C,G and T
bool is_good_kmer(const char* S, int64_t start, int64_t k){
    for(int64_t i = start; i < start + k; i++){
        if(S[i] != 'A' && S[i] != 'C' && S[i] != 'G' && S[i] != 'T') return false;
    }
    return true;
}

string get_rc(const string& S){
    string revc(S.rbegin(), S.rend());
    for(char& c : revc) c = rc_table[c];
    return revc;
}

template<typename output_stream_t>
void write_number_in_ascii(output_stream_t& out, int64_t x){
    string ascii = to_string(x);
    out.write(ascii.c_str(), ascii.size()); 
}