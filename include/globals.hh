#pragma once

#include <string>
#include <vector>

using namespace std;

void create_directory_if_does_not_exist(string path);

// https://stackoverflow.com/questions/18100097/portable-way-to-check-if-directory-exists-windows-linux-c
void check_dir_exists(string path);

// Returns filename of a new color file that has one color for each sequence
// Input format is either "fasta" or "fastq"
string generate_default_colorfile(string inputfile);

pair<string,string> split_all_seqs_at_non_ACGT(string inputfile, string inputfile_format, string colorfile);

vector<int64_t> read_colorfile(string filename);

std::string fix_alphabet(const std::string& input_file);

bool colex_compare(const string& S, const string& T);

// Split by whitespace
vector<string> split(string text);

// Split by delimiter
vector<string> split(string text, char delimiter);
vector<string> split(const char* text, char delimiter);

void reverse_complement_c_string(char* S, int64_t len);

bool files_are_equal(const std::string& p1, const std::string& p2);

template <typename T>
vector<T> parse_tokens(string S){
    vector<T> tokens;
    stringstream ss(S);
    T token;
    while(ss >> token) tokens.push_back(token);
    
    return tokens;
}

int64_t string_to_integer_safe(const string& S);

int64_t DNA_to_char_idx(char c); // ACGT -> 0123
char char_idx_to_DNA(int64_t i); // 0123 -> ACGT

void sigint_handler(int sig);

void sigabrt_handler(int sig);
