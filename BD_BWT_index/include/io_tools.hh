#ifndef TOOLS_HH
#define TOOLS_HH

#include <cstdlib>
#include <string>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

void copy_file(std::string file_from, std::string file_to);
void write_to_disk(std::string& text, std::string filepath);
char* read_from_disk_c_string(std::string filepath);
std::string read_from_disk(std::string filepath);


#endif
