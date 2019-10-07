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
#include <chrono>
#include <iomanip>

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


string fix_alphabet(string fastafile){
    write_log("Making all characters upper case and replacing non-{A,C,G,T} characters with 'A'");
    FASTA_reader fr(fastafile);
    string new_filename = temp_file_manager.get_temp_file_name("seqs-");
    ofstream out(new_filename);
    LL chars_replaced = 0;
    while(!fr.done()){
        string read = fr.get_next_query_stream().get_all();
        for(LL i = 0; i < read.size(); i++){
            char c = read[i];
            char c_new = toupper(c);
            if(c_new != 'A' && c_new != 'C' && c_new != 'G' && c_new != 'T'){
                c_new = 'A';
            }
            if(c_new != c){
                read[i] = c_new;
                chars_replaced++;
            }
        }
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
    FASTA_reader fr(fastafile);
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
    ifstream F(filename);
    if(!F.good()){
        cerr << "Error reading file: " << filename << endl;
        exit(1);
    }
}

// Also clears the file
void check_writable(string filename){
    ofstream F(filename, std::ofstream::out | std::ofstream::app);
    if(!F.good()){
        cerr << "Error writing to file: " << filename << endl;
        exit(1);
    }
}


template <typename T>
void write_to_file(string path, T& thing){
    ofstream out(path);
    if(!out.good()){
        cerr << "Error writing to " << path << endl;
        exit(1);
    }
    out << thing << endl;
}


template <typename T>
void read_from_file(string path, T& thing){
    ifstream input(path);
    if(!input.good()){
        cerr << "Error reading file: " << path << endl;
        exit(1);
    } else{
        input >> thing;
    }    
}

vector<string> get_all_lines(string infile){
    vector<string> lines;
    string line;
    ifstream in(infile);
    while(getline(in,line)){
        lines.push_back(line);
    }
    return lines;
}

vector<char> read_binary_file(string infile){
    std::ifstream file(infile, std::ios::binary | std::ios::ate);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

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
    std::ifstream f1(p1, std::ifstream::binary|std::ifstream::ate);
    std::ifstream f2(p2, std::ifstream::binary|std::ifstream::ate);

    if (f1.fail() || f2.fail()) {
        assert(false);
        return false; //file problem
    }

    if (f1.tellg() != f2.tellg()) {
      return false; //size mismatch
    }

    //seek back to beginning and use std::equal to compare contents
    f1.seekg(0, std::ifstream::beg);
    f2.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.rdbuf()));
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
