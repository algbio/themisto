#include "KallistoLite.hh"
#include "globals.hh"
#include <string>
#include <cstring>

using namespace std;

struct Config{
    LL k = -1;
    string fastafile;
    string colorfile;
    string index_dir;
    string temp_dir;
    bool load_boss = false;

    void check_valid(){
        if(fastafile != ""){
            check_readable(fastafile);
        }

        if(!load_boss){
            assert(k != -1);
        }

        if(colorfile != ""){
            check_readable(colorfile);
        }
        
        assert(index_dir != "");
        check_dir_exists(index_dir);

        assert(temp_dir != "");
        check_dir_exists(temp_dir);
    }
};

int main(int argc, char** argv){
    KallistoLite kl;
    if(argc == 1){
        cerr << "Options: " << endl;
        cerr << "  --load-boss (if given, loads a precomputed boss from the index directory)" << endl;
        cerr << "  --k [value of k] (required only if --load-boss is not given)" << endl;
        cerr << "  --fasta-file [filename] (always required)" << endl,
        cerr << "  --color-file [filename] (required only if you want to build the colors)" << endl;
        cerr << "  --index-dir [path] (always required, directory must exist before running)" << endl;
        cerr << "  --temp-dir [path] (always required, directory must exist before running)" << endl;
        cerr << "Usage examples:" << endl;
        cerr << "Build BOSS and colors:" << endl;
        cerr << "  ./build_index --k 31 --fasta-file references.fna --color-file colors.txt --index-dir index --temp-dir temp" << endl;
        cerr << "Build only the BOSS" << endl;
        cerr << "  ./build_index --k 31 --fasta-file references.fna --index-dir index --temp-dir temp" << endl;
        cerr << "Load a previously built BOSS from the index directory and compute the colors:" << endl;
        cerr << "  ./build_index --fasta-file references.fna --color-file colors.txt --index-dir index --temp-dir temp --load-boss" << endl;
        exit(1);
    }

    Config C;
    for(auto keyvalue : parse_args(argc, argv)){
        string option = keyvalue.first;
        vector<string> values = keyvalue.second;
        if(option == "--k"){
            assert(values.size() == 1);
            C.k = std::stoll(values[0]);
        } else if(option == "--fasta-file"){
            assert(values.size() == 1);
            C.fastafile = values[0];
        } else if(option == "--color-file"){
            assert(values.size() == 1);
            C.colorfile = values[0];
        } else if(option == "--index-dir"){
            assert(values.size() == 1);
            C.index_dir = values[0];
        } else if(option == "--temp-dir"){
            assert(values.size() == 1);
            C.temp_dir = values[0];
        } else if(option == "--load-boss"){
            C.load_boss = true;
        } else{
            cerr << "Error parsing command line arguments. Unkown option: " << option << endl;
            exit(1);
        }
    }

    C.check_valid();

    write_log("Starting");

    if(C.load_boss){
        write_log("Loading BOSS");
        kl.load_boss(C.index_dir + "/boss-");
    } else{
        write_log("Building BOSS");
        kl.construct_boss(C.fastafile, C.k);
        kl.save_boss(C.index_dir + "/boss-");
        write_log("Building BOSS finished (" + std::to_string(kl.boss.get_number_of_nodes()) + " nodes)");
    }

    if(C.colorfile != ""){
        write_log("Building colors");
        kl.construct_colors(C.fastafile, C.colorfile);
        kl.save_colors(C.index_dir + "/coloring-");
    }

    write_log("Finished");
}
