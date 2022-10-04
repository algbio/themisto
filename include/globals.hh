#pragma once

#include <string>
#include <vector>

using namespace std;

void create_directory_if_does_not_exist(string path);

// https://stackoverflow.com/questions/18100097/portable-way-to-check-if-directory-exists-windows-linux-c
void check_dir_exists(string path);

// Returns filename of a new color file that has one color for each sequence
// Input format is either "fasta" or "fastq"
string generate_default_colorfile(string inputfile, string file_format);

pair<string,string> split_all_seqs_at_non_ACGT(string inputfile, string inputfile_format, string colorfile);

vector<int64_t> read_colorfile(string filename);

std::string fix_alphabet(const std::string& input_file);