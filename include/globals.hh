#pragma once


#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cmath>
#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include <string>
#include <map>

typedef int64_t LL;
static const char read_separator = '$';
#ifndef KMER_MAX_LENGTH
#define KMER_MAX_LENGTH 64
#endif

#include "TempFileManager.hh"
#include <signal.h>
#include "input_reading.hh"
#include "throwing_streams.hh"
#include <chrono>
#include <iomanip>
#include <random>
#include <filesystem>

using namespace std;
using namespace std::chrono;

// Returns a reference to the singleton temp file manager
Temp_File_Manager& get_temp_file_manager();

long long cur_time_millis();
double seconds_since_program_start();
string getTimeString();
void enable_logging();
void disable_logging();
void write_log(string message);
map<string,vector<string> > parse_args(int argc, char** argv);
string figure_out_file_format(string filename);
char fix_char(char c);

// Returns number of chars replaced
LL fix_alphabet_of_string(string& S);

// Makes a copy of the file and replaces bad characters. Returns the new filename
// The new file is in fasta format
std::string fix_alphabet(const std::string& input_file, const std::size_t bufsiz, const int mode);

// We need this function because the standard library stoll function accepts all kinds of crap,
// such as "123aasfjhk" and "1 2 3 4" as a number. This function check that the string is a valid
// number and returns that number, or throws an error otherwise.
LL string_to_integer_safe(const string& S);

vector<LL> read_colorfile(string filename);

// Returns new inputfile and new colorfile
pair<string,string> split_all_seqs_at_non_ACGT(string inputfile, string inputfile_format, string colorfile);

void sigint_handler(int sig);
void sigabrt_handler(int sig);

vector<string> get_first_and_last_kmers(string fastafile, LL k);

string get_rc(string S);

// true if S is colexicographically-smaller than T
bool colex_compare(const string& S, const string& T);

bool colex_compare_cstrings(const char* x, const char* y);
bool lex_compare(const string& S, const string& T);
bool lex_compare_cstrings(const char* x, const char* y);

// Split by whitespace
vector<string> split(string text);

// Split by delimiter
vector<string> split(string text, char delimiter);
vector<string> split(const char* text, char delimiter);

void create_directory_if_does_not_exist(string path);
void check_dir_exists(string path);
void check_readable(string filename);

// Also clears the file
void check_writable(string filename);

vector<string> get_all_lines(string infile);
vector<char> read_binary_file(string infile);
bool files_are_equal(const std::string& p1, const std::string& p2);
void check_true(bool condition, string error_message);

// Returns filename of a new color file that has one color for each sequence
// Input format is either "fasta" or "fastq"
string generate_default_colorfile(string inputfile, string file_format);

template <typename T>
vector<T> parse_tokens(string S){
    vector<T> tokens;
    stringstream ss(S);
    T token;
    while(ss >> token) tokens.push_back(token);
    
    return tokens;
}

// Returns filename 
string string_to_temp_file(const string& S);

template <typename T>
void write_to_file(string path, T& thing){
    throwing_ofstream out(path);
    out << thing << endl;
}

template <typename T>
void read_from_file(string path, T& thing){
    throwing_ifstream input(path);
    input >> thing;
}

class Progress_printer{

    public:

    LL n_jobs;
    LL processed;
    LL total_prints;
    LL next_print;

    Progress_printer(LL n_jobs, LL total_prints) : n_jobs(n_jobs), processed(0), total_prints(total_prints), next_print(0) {}

    void job_done(){
        if(next_print == processed){
            LL progress_percent = round(100 * ((double)processed / n_jobs));
            write_log("Progress: " + to_string(progress_percent) + "%");
            next_print += n_jobs / total_prints;
        }
        processed++;
    }

};

