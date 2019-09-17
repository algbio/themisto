#include "io_tools.hh"
#include <string>
#include <cassert>
#include <ctime>
#include <chrono>
#include <fstream>
#include <cstdio>
#include <cerrno>

using namespace std;

void write_to_disk(string& text, string filepath){
    ofstream stream(filepath);
    stream << text;
    stream.flush();
}

char* read_from_disk_c_string(std::string filename)
{
  std::FILE *fp = std::fopen(filename.c_str(), "rb");
  if (fp)
  {
    std::fseek(fp, 0, SEEK_END);
    int64_t size = std::ftell(fp);
    char* contents = (char*)malloc(size + 1);
    std::rewind(fp);
    int unused = std::fread(&contents[0], 1, size, fp);
    (void) unused; // TODO: do something with this value
    std::fclose(fp);
    return contents;
  }
  throw(errno);
}

std::string read_from_disk(std::string filename)
{
  std::FILE *fp = std::fopen(filename.c_str(), "rb");
  if (fp)
  {
    std::string contents;
    std::fseek(fp, 0, SEEK_END);
    contents.resize(std::ftell(fp));
    std::rewind(fp);
    int unused = std::fread(&contents[0], 1, contents.size(), fp);
    (void) unused; // TODO: do something with this value
    std::fclose(fp);
    return contents;
  }
  throw(errno);
}

void copy_file(std::string file_from, std::string file_to){
    ifstream from(file_from);
    ofstream to(file_to);
    string line;
    while(getline(from, line)) 
        to << line << "\n";
}