#include "Themisto.hh"
#include "input_reading.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"

using namespace std;

struct Config{
    vector<string> query_files;
    vector<string> outfiles;
    string index_dir;
    string temp_dir;
    
    bool gzipped_output = false;
    bool reverse_complements = false;
    bool sort_output = false;
    LL n_threads = 1;

    void check_valid(){
        check_true(index_dir != "", "Index directory not set");
        check_dir_exists(index_dir);

        check_true(temp_dir != "", "Temp directory not set");
        check_dir_exists(temp_dir);
    }
};

vector<string> read_lines(string filename){
    check_readable(filename);
    vector<string> lines;
    throwing_ifstream in(filename);
    string line;
    while(in.getline(line)){
        lines.push_back(line);
    }
    return lines;
}

int main2(int argc, char** argv){
    if(argc != 3){
        cerr << "Usage: " << argv[0] << " index_directory temp_directory" << endl;
        cerr << "Prints to stdout one line for each node. The line contains the colors of the node as space-separated integers" << endl;
        return 1;
    }
    Config C;
    C.index_dir = argv[1];
    C.temp_dir = argv[2];
    C.check_valid();

    write_log("Starting");
    temp_file_manager.set_dir(C.temp_dir);

    write_log("Loading the index");    
    Themisto themisto;
    themisto.load_boss(C.index_dir + "/boss-");
    themisto.load_colors(C.index_dir + "/coloring-");

    write_log("Dumping colors");    

    vector<LL> colors(themisto.coloring.n_colors); // List of colors
    for(LL i = 0; i < themisto.boss.number_of_nodes(); i++){
        LL colorset_size = themisto.coloring.get_colorset_to_buffer(i, themisto.boss, colors);
        for(LL j = 0; j < colorset_size; j++) cout << colors[j] << " ";
        cout << "\n";
    }

    write_log("Finished");

    return 0;
}

int main(int argc, char** argv){
    write_log("Themisto-" + std::string(THEMISTO_BUILD_VERSION));
    try{
        return main2(argc, argv);
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    }
}
