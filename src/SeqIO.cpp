#include "SeqIO.hh"
#include <algorithm>

namespace SeqIO{

using namespace SeqIO;

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

} // namespace SeqIO
