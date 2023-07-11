#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "commands.hh"
#include "sbwt/globals.hh"
#include "version.h"

using namespace std;

static vector<string> commands = {"build", "pseudoalign", "extract-unitigs", "dump-color-matrix", "stats"};

void print_help(int argc, char** argv){
    (void) argc; // Unused parameter
    cerr << "Available commands: " << endl;
    for(string S : commands) cerr << "   " << argv[0] << " " << S << endl;
    cerr << "Running a command without arguments prints the usage instructions for the command." << endl;
}

int main(int argc, char** argv){

    bool silent_in_argv = false;
    for(int64_t i = 1; i < argc; i++)
        if(string(argv[i]) == string("--silent"))
            silent_in_argv = true;

    if(!silent_in_argv){
        sbwt::write_log("Themisto-" + std::string(THEMISTO_BUILD_VERSION), sbwt::LogLevel::MAJOR);
        sbwt::write_log("Maximum k-mer length (size of the de Bruijn graph node labels): " + std::to_string(MAX_KMER_LENGTH), sbwt::LogLevel::MAJOR);
    }

    if(argc == 1){
        print_help(argc, argv);
        return 0;
    }

    string command = argv[1];
    if(command == "--help" || command == "-h"){
        print_help(argc, argv);
        return 0;
    }

    // Drop the first element of argv
    for(int64_t i = 1; i < argc; i++) argv[i-1] = argv[i];
    argc--;

    try{
        if(command == "build") return build_index_main(argc, argv);
        else if(command == "pseudoalign") return pseudoalign_main(argc, argv);
        else if(command == "extract-unitigs") return extract_unitigs_main(argc, argv);
        else if(command == "stats") return stats_main(argc, argv);
        else if(command == "dump-color-matrix") return dump_color_matrix_main(argc, argv); // Undocumented developer feature
        else if(command == "color-set-diagnostics") return color_set_diagnostics_main(argc, argv); // Undocumented developer feature
        else if(command == "make-d-equal-1") return make_d_equal_1_main(argc, argv); // Undocumented developer feature
        else if(command == "dump-distinct-color-sets-to-binary") return dump_distinct_color_sets_to_binary(argc, argv); // Undocumented developer feature
        else{
            throw std::runtime_error("Invalid command: " + command);
            return 1;
        }
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    } catch(const std::exception& e){
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

}