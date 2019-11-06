#include "KallistoLite.hh"
#include "globals.hh"
#include <string>
#include <cstring>

using namespace std;

struct Config{
    LL k = -1;
    LL n_threads = 1;
    string fastafile;
    string colorfile;
    string index_dir;
    string temp_dir;
    bool load_boss = false;
    LL memory_megas = 1000;
    bool auto_colors = false;

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

    string to_string(){
        stringstream ss;
        ss << "Input fasta-file = " << fastafile << "\n";
        if(colorfile != "") ss << "Color name file = " << colorfile << "\n";
        ss << "Index directory = " << index_dir << "\n";
        ss << "Temporary directory = " << temp_dir << "\n";
        ss << "k = " << k << "\n";
        ss << "Number of threads = " << n_threads << "\n";
        ss << "Memory megabytes = " << memory_megas << "\n";
        ss << "Automatic colors = " << (auto_colors ? "true" : "false") << "\n";
        ss << "Load BOSS = " << (load_boss ? "true" : "false"); // Last has no endline
        return ss.str();
    }
};

int main2(int argc, char** argv){
    KallistoLite kl;
    if(argc == 1){
        cerr << "Options: " << endl;
        cerr << "  --load-boss (if given, loads a precomputed boss from the index directory)" << endl;
        cerr << "  --k [value of k] (required only if --load-boss is not given)" << endl;
        cerr << "  --fasta-file [filename] (always required)" << endl,
        cerr << "  --color-file [filename] (required only if you want to build the colors)" << endl;
        cerr << "  --auto-colors (instead of color file let the program automatically give colors integer names)" << endl;
        cerr << "  --index-dir [path] (always required, directory must exist before running)" << endl;
        cerr << "  --temp-dir [path] (always required, directory must exist before running)" << endl;
        cerr << "  --mem-megas [number] (Number of megabytes allowed for external memory algorithms. Default: 1000)" << endl;
        cerr << "  --disk-megas [number] (megabytes of disk space allocated for stxll library. Default: 10000)" << endl;
        cerr << "  --n-threads [number] (number of parallel threads to use. Default: 1)" << endl;
        cerr << "Usage examples:" << endl;
        cerr << "Build BOSS and colors:" << endl;
        cerr << "  ./build_index --k 31 --mem-megas 10000 --fasta-file references.fna --color-file colors.txt --index-dir index --temp-dir temp" << endl;
        cerr << "Build only the BOSS" << endl;
        cerr << "  ./build_index --k 31 --mem-megas 10000 --fasta-file references.fna --index-dir index --temp-dir temp" << endl;
        cerr << "Load a previously built BOSS from the index directory and compute the colors:" << endl;
        cerr << "  ./build_index --mem-megas 10000 --fasta-file references.fna --color-file colors.txt --index-dir index --temp-dir temp --load-boss" << endl;
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
        } else if(option == "--n-threads"){
            assert(values.size() == 1);
            C.n_threads = std::stoll(values[0]);
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
            assert(values.size() == 0);
            C.load_boss = true;
        } else if(option == "--mem-megas"){
            assert(values.size() == 1);
            C.memory_megas = std::stoll(values[0]);
        } else if(option == "--auto-colors"){
            assert(values.size() == 0);
            C.auto_colors = true;
        } else{
            cerr << "Error parsing command line arguments. Unkown option: " << option << endl;
            exit(1);
        }
    }

    C.check_valid();
    temp_file_manager.set_dir(C.temp_dir);

    cerr << C.to_string() << endl;
    write_log("Starting");
    C.fastafile = fix_FASTA_alphabet(C.fastafile);

    if(C.load_boss){
        write_log("Loading BOSS");
        kl.load_boss(C.index_dir + "/boss-");
    } else{
        write_log("Building BOSS");
        kl.construct_boss(C.fastafile, C.k, C.memory_megas * 1e6, C.n_threads);
        kl.save_boss(C.index_dir + "/boss-");
        write_log("Building BOSS finished (" + std::to_string(kl.boss.get_number_of_nodes()) + " nodes)");
    }

    if(C.colorfile != "" || C.auto_colors){
        write_log("Building colors");
        kl.construct_colors(C.fastafile, C.auto_colors ? "" : C.colorfile, C.memory_megas * 1e6, C.n_threads);
        kl.save_colors(C.index_dir + "/coloring-");
    }

    write_log("Finished");

    return 0;
}

int main(int argc, char** argv){
    try{
        return main2(argc, argv);
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    }
}