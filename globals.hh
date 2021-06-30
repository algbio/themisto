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
#include "TempFileManager.hh"
#include <signal.h>
#include "input_reading.hh"
#include "throwing_streams.hh"
#include <chrono>
#include <iomanip>
#include <random>

#ifndef KMER_MAX_LENGTH
#define KMER_MAX_LENGTH 64
#endif

using namespace std;
using namespace std::chrono;

typedef int64_t LL;
const char read_separator = '$';
Temp_File_Manager temp_file_manager;

long long cur_time_millis(){
	return (std::chrono::duration_cast< milliseconds >(system_clock::now().time_since_epoch())).count();
}

static long long int program_start_millis = cur_time_millis();

double seconds_since_program_start(){
	return (cur_time_millis() - program_start_millis) / 1000.0;
}


string getTimeString(){
    std::time_t result = std::time(NULL);
    string time = std::asctime(std::localtime(&result));
    return time.substr(0,time.size() - 1); // Trim the trailing newline
}

static bool logging_enabled = true;

void enable_logging(){
    logging_enabled = true;
}

void disable_logging(){
    logging_enabled = false;
}

std::mutex write_log_mutex;

void write_log(string message){
    std::lock_guard<std::mutex> lock(write_log_mutex);
    if(logging_enabled){
        std::streamsize default_precision = std::cout.precision();

        std::cerr << 
        std::setprecision(4) << std::fixed <<
        seconds_since_program_start() <<
        std::setprecision(default_precision) << 
        " " << getTimeString() << " " << message << std::endl;
    }
}

map<string,vector<string> > parse_args(int argc, char** argv){
    // Options are argumenta that start with "--". All non-options
    // that come after an option are parameters for that option
    map<string,vector<string> > M; // Option -> list of parameters
    string current_option = "";
    for(LL i = 1; i < argc; i++){
        string S = argv[i];
        if(S.size() >= 2  && S.substr(0,2) == "--"){
            current_option = S;
            M[current_option].resize(0); // Add empty vector for this option.
        } else{
            if(current_option == ""){
                cerr << "Error parsing command line parameters" << endl;
                exit(1);
            }
            M[current_option].push_back(S);
        }
    }
    return M;
}

string figure_out_file_format(string filename){
    for(LL i = filename.size()-1; i >= 0; i--){
        if(filename[i] == '.'){
            string end = filename.substr(i);
            
            if(end == ".fasta") return "fasta";
            if(end == ".fna") return "fasta";
            if(end == ".ffn") return "fasta";
            if(end == ".faa") return "fasta";
            if(end == ".frn") return "fasta";
            if(end == ".fa") return "fasta";

            if(end == ".fastq") return "fastq";
            if(end == ".fq") return "fastq";

            if(end == ".gz") return "gzip";

            throw(runtime_error("Unknown file format: " + filename));
        }
    }
    throw(runtime_error("Unknown file format: " + filename));
    return "unknown";
}

static constexpr char R_conv_tbl[] = { 'A', 'G' };
static constexpr char Y_conv_tbl[] = { 'C', 'T' };
static constexpr char K_conv_tbl[] = { 'G', 'T' };
static constexpr char M_conv_tbl[] = { 'A', 'C' };
static constexpr char S_conv_tbl[] = { 'C', 'G' };
static constexpr char W_conv_tbl[] = { 'A', 'T' };
static constexpr char B_conv_tbl[] = { 'C', 'G', 'T' };
static constexpr char D_conv_tbl[] = { 'A', 'G', 'T' };
static constexpr char H_conv_tbl[] = { 'A', 'C', 'T' };
static constexpr char V_conv_tbl[] = { 'A', 'C', 'G' };
static constexpr char N_conv_tbl[] = { 'A', 'C', 'G', 'T' };

char fix_char(char c){
	c = toupper(c);
    int rd = std::rand();
    
	switch (c) {
	case 'A':
		return c;
	case 'C':
		return c;
	case 'G':
		return c;
	case 'T':
		return c;
	case 'U':
		return 'T';
	case 'R':
		return R_conv_tbl[rd % 2];
	case 'Y':
		return Y_conv_tbl[rd % 2];
	case 'K':
		return K_conv_tbl[rd % 2];
	case 'M':
		return M_conv_tbl[rd % 2];
	case 'S':
		return S_conv_tbl[rd % 2];
	case 'W':
		return W_conv_tbl[rd % 2];
	case 'B':
		return B_conv_tbl[rd % 3];
	case 'D':
		return D_conv_tbl[rd % 3];
	case 'H':
		return H_conv_tbl[rd % 3];
	case 'V':
		return V_conv_tbl[rd % 3];
	default:
		return N_conv_tbl[rd % 4];
	}
}

// Returns number of chars replaced
LL fix_alphabet_of_string(string& S){
    LL chars_replaced = 0;
    for(LL i = 0; i < S.size(); i++){
        char c = S[i];
        char c_new = fix_char(c);
        if(c_new != c){
            S[i] = c_new;
            chars_replaced++;
        }
    }
    return chars_replaced;
}


// Makes a copy of the file and replaces a bad characters. Returns the new filename
// The new file is in fasta format
string fix_alphabet(Sequence_Reader& sr){
    write_log("Making all characters upper case and replacing non-{A,C,G,T} characters with random characeters from {A,C,G,T}");
    //Sequence_Reader fr(fastafile, FASTA_MODE);
    string new_filename = temp_file_manager.get_temp_file_name("seqs-");
    throwing_ofstream out(new_filename);
    LL chars_replaced = 0;
    while(!sr.done()){
        string read = sr.get_next_query_stream().get_all();
        chars_replaced += fix_alphabet_of_string(read);
        out << ">\n" << read << "\n";
    }
    write_log("Replaced " + to_string(chars_replaced) + " characters");
    return new_filename;
}

void sigint_handler(int sig) {
    cerr << "caught signal: " << sig << endl;
    cerr << "Cleaning up temporary files" << endl;
    temp_file_manager.clean_up();
    exit(1);
}

void sigabrt_handler(int sig) {
    cerr << "caught signal: " << sig << endl;
    cerr << "Cleaning up temporary files" << endl;
    temp_file_manager.clean_up();
    cerr << "Aborting" << endl;
    exit(1);
}

auto sigint_register_return_value = signal(SIGINT, sigint_handler); // Set the SIGINT handler
auto sigabrt_register_return_value = signal(SIGABRT, sigabrt_handler); // Set the SIGABRT handler

void set_temp_dir(string temp_dir){
    temp_file_manager.set_dir(temp_dir);
}

string get_temp_file_name(string prefix){
    return temp_file_manager.get_temp_file_name(prefix);
}

vector<string> get_first_and_last_kmers(string fastafile, LL k){
    // todo: this is pretty expensive because this has to read the whole reference data
    Sequence_Reader fr(fastafile, FASTA_MODE);
    vector<string> result;
    while(!fr.done()){
        string ref = fr.get_next_query_stream().get_all();
        if(ref.size() >= k){
            result.push_back(ref.substr(0,k));
            result.push_back(ref.substr(ref.size()-k,k));
        }
    }
    return result;
}

string get_rc(string S){
    std::reverse(S.begin(), S.end());
    for(char& c : S){
        if(c == 'A') c = 'T';
        else if(c == 'C') c = 'G';
        else if(c == 'G') c = 'C';
        else if(c == 'T') c = 'A';
    }
    return S;
}

// true if S is colexicographically-smaller than T
bool colex_compare(const string& S, const string& T){
    LL i = 0;
    while(true){
        if(i == S.size() || i == T.size()){
            // One of the strings is a suffix of the other. Return the shorter.
            if(S.size() < T.size()) return true;
            else return false;
        }
        if(S[S.size()-1-i] < T[T.size()-1-i]) return true;
        if(S[S.size()-1-i] > T[T.size()-1-i]) return false;
        i++;
    }
}

bool colex_compare_cstrings(const char* x, const char* y){
    LL nx = strlen(x);
    LL ny = strlen(y);
    for(LL i = 0; i < min(nx,ny); i++){
        if(x[nx-1-i] < y[ny-1-i]) return true;
        if(x[nx-1-i] > y[ny-1-i]) return false;
    }

    // All no mismatches -> the shorter string is smaller
    return nx < ny;
};

bool lex_compare(const string& S, const string& T){
    return S < T;
};

bool lex_compare_cstrings(const char* x, const char* y){
    return strcmp(x,y) < 0;
};


template <typename T>
vector<T> parse_tokens(string S){
    vector<T> tokens;
    stringstream ss(S);
    T token;
    while(ss >> token) tokens.push_back(token);
    
    return tokens;
}

// Split by whitespace
vector<string> split(string text){
    std::istringstream iss(text);
    std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                 std::istream_iterator<std::string>());
    return results;
}

// Split by delimiter
vector<string> split(string text, char delimiter){
    assert(text.size() != 0); // If called with empty string we probably have a bug
    vector<LL> I; // Delimiter indices
    I.push_back(-1);
    for(LL i = 0; i < text.size(); i++){
        if(text[i] == delimiter){
            I.push_back(i);
        }
    }
    I.push_back(text.size());
    vector<string> tokens;
    for(LL i = 0; i < I.size()-1; i++){
        LL len = I[i+1] - I[i] + 1 - 2;
        tokens.push_back(text.substr(I[i]+1, len));
    }
    
    return tokens;
}

vector<string> split(const char* text, char delimiter){
    return split(string(text), delimiter);
}


// https://stackoverflow.com/questions/18100097/portable-way-to-check-if-directory-exists-windows-linux-c
void check_dir_exists(string path){
    struct stat info;    
    if( stat( path.c_str(), &info ) != 0 ){
        cerr << "Error: can not access directory " << path << endl;
        exit(1);
    }
    else if( info.st_mode & S_IFDIR ){
        // All good
    }    
    else{
        cerr << "Error: is not a directory: " << path << endl;
        exit(1);
    }
}

void check_readable(string filename){
    throwing_ifstream F(filename); // Throws on failure
}

// Also clears the file
void check_writable(string filename){
    throwing_ofstream F(filename, std::ofstream::out | std::ofstream::app); // Throws on failure
}


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

vector<string> get_all_lines(string infile){
    vector<string> lines;
    string line;
    throwing_ifstream in(infile);
    while(in.getline(line)){
        lines.push_back(line);
    }
    return lines;
}

vector<char> read_binary_file(string infile){
    throwing_ifstream file(infile, std::ios::binary | std::ios::ate);
    std::streamsize size = file.stream.tellg();
    file.stream.seekg(0, std::ios::beg);

    std::vector<char> buffer(size);
    if (file.read(buffer.data(), size)){
        return buffer;
    } else{
        cerr << "Error reading file: " << infile << endl;
        assert(false);
    }
}

bool files_are_equal(const std::string& p1, const std::string& p2) {
  //https://stackoverflow.com/questions/6163611/compare-two-files/6163627
    throwing_ifstream f1(p1, std::ifstream::binary|std::ifstream::ate);
    throwing_ifstream f2(p2, std::ifstream::binary|std::ifstream::ate);

    if (f1.stream.tellg() != f2.stream.tellg()) {
      return false; //size mismatch
    }

    //seek back to beginning and use std::equal to compare contents
    f1.stream.seekg(0, std::ifstream::beg);
    f2.stream.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(f1.stream.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.stream.rdbuf()));
}

void check_true(bool condition, string error_message){
    if(!condition){
        throw std::runtime_error(error_message);
    }
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

