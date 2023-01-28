#include "SeqIO.hh"
#include <algorithm>

namespace SeqIO{

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


char get_rc(char c){
    return rc_table[(unsigned char)c];
}

string get_rc(const string& S){
    string T = S;
    std::reverse(T.begin(), T.end());
    for(char& c : T) c = get_rc(c);
    return T;
}

const vector<string> fasta_suffixes = {".fasta", ".fna", ".ffn", ".faa", ".frn", ".fa"};
const vector<string> fastq_suffixes = {".fastq", ".fq"};

FileFormat figure_out_file_format(string filename){
    Format fasta_or_fastq;
    bool gzipped = false;
    string extension;

    string gzip_suffix = "";
    if(filename.size() >= 3 && filename.substr(filename.size()-3) == ".gz"){
        filename = filename.substr(0, filename.size()-3); // Drop .gz
        gzipped = true;
    }

    for(int64_t i = (int64_t)filename.size()-1; i >= 0; i--){
        if(filename[i] == '.'){
            string ending = filename.substr(i);
            
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

            throw(runtime_error("Unknown file format: " + filename + (gzipped ? ".gz" : "")));
        }
    }
    throw(runtime_error("Unknown file format: " + filename + (gzipped ? ".gz" : "")));
}

vector<uint8_t> interpret_fastq_quality_scores(const char* qual, int64_t len){
    vector<uint8_t> v(len);
    for(int64_t i = 0; i < len; i++){
        v[i] = qual[i] - 0x21;
    }
    return v;
}

} // namespace SeqIO