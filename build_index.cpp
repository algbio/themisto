#include "Themisto.hh"
#include "globals.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"

using namespace std;

struct Config{
    LL k = -1;
    LL n_threads = 1;
    string inputfile;
    string colorfile;
    string index_dir;
    string temp_dir;
    string input_format;
    bool load_boss = false;
    LL memory_megas = 1000;
    bool auto_colors = false;
    bool del_non_ACGT = false;
    LL pp_buf_siz = 1024*4;
    LL colorset_sampling_distance = 1;
    
    void check_valid(){
        check_true(inputfile != "", "Input file not set");

        check_readable(inputfile);
        check_true(input_format != "", "Problem detecting input format");

        if(!load_boss){
            check_true(k != -1, "Parameter k not set");
            check_true(k <= KMER_MAX_LENGTH, "Maximum allowed k is " + std::to_string(KMER_MAX_LENGTH) + ". To increase the limit, recompile by first running cmake with the option `-DMAX_KMER_LENGTH=n`, where n is a number up to 255, and then running `make` again."); // 255 is max because of KMC
        }

        if(colorfile != ""){
            check_readable(colorfile);
        }
        
        check_true(index_dir != "", "Index directory not set");
        check_dir_exists(index_dir);

        check_true(temp_dir != "", "Temp directory not set");
        check_dir_exists(temp_dir);

        check_true(memory_megas > 0, "Memory budget must be positive");
        check_true(colorset_sampling_distance >= 1, "Colorset sampling distance must be positive");
        check_true(pp_buf_siz > 0, "Preprocessing buffer size must be positive");

    }

    string to_string(){
        stringstream ss;
        ss << "Input file = " << inputfile << "\n";
        ss << "Input format = " << input_format << "\n";
        if(colorfile != "") ss << "Color name file = " << colorfile << "\n";
        ss << "Index directory = " << index_dir << "\n";
        ss << "Temporary directory = " << temp_dir << "\n";
        ss << "k = " << k << "\n";
        ss << "Number of threads = " << n_threads << "\n";
        ss << "Memory megabytes = " << memory_megas << "\n";
        ss << "Automatic colors = " << (auto_colors ? "true" : "false") << "\n";
        ss << "Load BOSS = " << (load_boss ? "true" : "false") << "\n";
        ss << "Handing of non-ACGT characters = " << (del_non_ACGT ? "Delete" : "Randomize") << "\n";
        ss << "Preprocessing buffer size = " << pp_buf_siz; // Last has no endline
        return ss.str();
    }
};

// Returns filename of a new color file that has one color for each sequence
// Input format is either "fasta" or "fastq"
string generate_default_colorfile(string inputfile, string file_format){
    string colorfile = temp_file_manager.get_temp_file_name("");
    throwing_ofstream out(colorfile);
    Sequence_Reader fr(inputfile, file_format == "fasta" ? FASTA_MODE : FASTQ_MODE);
    LL seq_id = 0;
    while(!fr.done()){
        fr.get_next_query_stream().get_all();
        out << seq_id << "\n";
        seq_id++;
    }
    return colorfile;
}

int main2(int argc, char** argv){

    // Legacy support: transform old option format --k to -k
    string legacy_support_fix = "-k";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--k") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "Builds an index consisting of compact de Bruijn graph using the BOSS data structure and color information. The input is a set of reference sequences in a single file in fasta or fastq format, and a colorfile, which is a plain text file containing the colors of the reference sequences in the same order as they appear in the reference sequence file, one line per sequence. The names are given as ASCII strings, but they should not contain whitespace characters. If there are characters outside of the DNA alphabet ACGT in the input sequences, those are replaced with random characters from the DNA alphabet.");

    options.add_options()
        ("load-boss", "If given, loads a precomputed BOSS from the index directory", cxxopts::value<bool>()->default_value("false"))
        ("k,node-length", "The k of the k-mers. Required only if --load-boss is not given", cxxopts::value<LL>())
        ("i,input-file", "The input sequences in FASTA or FASTQ format. The format is inferred from the file extension. Recognized file extensions for fasta are: .fasta, .fna, .ffn, .faa and .frn . Recognized extensions for fastq are: .fastq and .fq . If the file ends with .gz, it is uncompressed into a temporary directory and the temporary file is deleted after use.", cxxopts::value<string>())
        ("c,color-file", "One color per sequence in the fasta file, one color per line. If not given, colors are not built, unless --auto-colors is given.", cxxopts::value<string>()->default_value(""))
        ("auto-colors", "Instead of a color file, number the sequences with integers 0,1,2,... in the same order as in the sequence file.)", cxxopts::value<bool>()->default_value("false"))
        ("o,index-dir", "Directory where the index will be built. Always required, directory must exist before running.", cxxopts::value<string>())
        ("d,colorset-pointer-tradeoff", "This option controls a time-space tradeoff for storing and querying color sets. If given a value d, we store color set pointers only for every d nodes on every unitig. The higher the value of d, the smaller then index, but the slower the queries. The savings might be significant if the number of distinct color sets is small and the graph is large and has long unitigs.", cxxopts::value<LL>()->default_value("1"))
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
        ("m,mem-megas", "Number of megabytes allowed for external memory algorithms. Default: 1000", cxxopts::value<LL>()->default_value("1000"))
        ("t, n-threads", "Number of parallel exectuion threads. Default: 1", cxxopts::value<LL>()->default_value("1"))
        ("delete-non-ACGT", "Delete k-mers that have a letter outside of the DNA alphabet ACGT. If this option is not given, the non-ACGT letters are replaced with random nucleotides.", cxxopts::value<bool>()->default_value("false"))
        ("pp-buf-siz", "Size of preprocessing buffer (in bytes) for fixing alphabet", cxxopts::value<LL>()->default_value("4096"))
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        cerr << "Usage examples:" << endl;
        cerr << "Build BOSS and colors:" << endl;
        cerr << "  ./build_index -k 31 --mem-megas 10000 --input-file references.fna --color-file colors.txt --index-dir index --temp-dir temp" << endl;
        cerr << "Build only the BOSS" << endl;
        cerr << "  ./build_index -k 31 --mem-megas 10000 --input-file references.fna --index-dir index --temp-dir temp" << endl;
        cerr << "Load a previously built BOSS from the index directory and compute the colors:" << endl;
        cerr << "  ./build_index --mem-megas 10000 --input-file references.fna --color-file colors.txt --index-dir index --temp-dir temp --load-boss" << endl;
        exit(1);
    }    

    Config C;
    Themisto themisto;
    C.k = opts["k"].as<LL>();
    C.inputfile = opts["input-file"].as<string>();
    C.input_format = figure_out_file_format(C.inputfile);
    C.n_threads = opts["n-threads"].as<LL>();
    C.colorfile = opts["color-file"].as<string>();
    C.index_dir = opts["index-dir"].as<string>();
    C.temp_dir = opts["temp-dir"].as<string>();
    C.load_boss = opts["load-boss"].as<bool>();
    C.memory_megas = opts["mem-megas"].as<LL>();
    C.auto_colors = opts["auto-colors"].as<bool>();
    C.pp_buf_siz = opts["pp-buf-siz"].as<LL>();
    C.colorset_sampling_distance = opts["colorset-pointer-tradeoff"].as<LL>();
    C.del_non_ACGT = opts["delete-non-ACGT"].as<bool>();

    create_directory_if_does_not_exist(C.index_dir);
    create_directory_if_does_not_exist(C.temp_dir);

    C.check_valid();
    temp_file_manager.set_dir(C.temp_dir);

    cerr << C.to_string() << endl;
    write_log("Starting");

    if(C.input_format == "gzip"){
        write_log("Decompressing the input file");
        string new_name = temp_file_manager.get_temp_file_name("input");
        check_true(gz_decompress(C.inputfile, new_name) == Z_OK, "Problem with zlib decompression");
        C.input_format = figure_out_file_format(C.inputfile.substr(0,C.inputfile.size() - 3));
        C.inputfile = new_name;
    }

    if(C.colorfile == ""){
        // Automatic colors
        C.colorfile = generate_default_colorfile(C.inputfile, C.input_format);
    }

    // Deal with non-ACGT characters

    Sequence_Reader sr(C.inputfile, C.input_format == "fasta" ? FASTA_MODE : FASTQ_MODE);
    if(C.del_non_ACGT){
        std::tie(C.inputfile, C.colorfile) = split_all_seqs_at_non_ACGT(C.inputfile, C.input_format, C.colorfile); // Turns the file into fasta format also
        C.input_format = "fasta"; // split_all_seqs_at_non_ACGT returns a fasta file
    } else {
        C.inputfile = fix_alphabet(C.inputfile, C.pp_buf_siz, C.input_format == "fasta" ? FASTA_MODE : FASTQ_MODE); // Turns the file into fasta format also
        C.input_format = "fasta"; // fix_alphabet returns a fasta file
    }
    
    if(C.load_boss){
        write_log("Loading BOSS");
        themisto.load_boss(C.index_dir + "/boss-");
    } else{
        write_log("Building BOSS");
        themisto.construct_boss(C.inputfile, C.k, C.memory_megas * 1e6, C.n_threads, false);
        themisto.save_boss(C.index_dir + "/boss-");
        write_log("Building BOSS finished (" + std::to_string(themisto.boss.number_of_nodes()) + " nodes)");
    }

    if(C.colorfile != "" || C.auto_colors){
        write_log("Building colors");
        themisto.construct_colors(C.inputfile, C.auto_colors ? "" : C.colorfile, C.memory_megas * 1e6, C.n_threads, C.colorset_sampling_distance);
        themisto.save_colors(C.index_dir + "/coloring-");
    }

    write_log("Finished");

    return 0;
}

int main(int argc, char** argv){
    write_log("Themisto-" + std::string(THEMISTO_BUILD_VERSION));
    write_log("Maximum k-mer length (size of the de Bruijn graph node labels): " + std::to_string(KMER_MAX_LENGTH-1));
    try{
        return main2(argc, argv);
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        return 1;
    } catch(const std::exception& e){
        std::cerr << "Error: " << e.what() << '\n';
    }
}
