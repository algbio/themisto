#include "Themisto.hh"
#include "globals.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"

using namespace std;

struct Build_Config{
    LL k = -1;
    LL n_threads = 1;
    string inputfile;
    string colorfile;
    string index_dbg_file;
    string index_color_file;
    string temp_dir;
    string input_format;
    bool load_dbg = false;
    LL memory_megas = 1000;
    bool no_colors = false;
    bool del_non_ACGT = false;
    LL colorset_sampling_distance = 1;
    
    void check_valid(){
        check_true(inputfile != "", "Input file not set");

        check_readable(inputfile);
        check_true(input_format != "", "Problem detecting input format");

        if(!load_dbg){
            check_true(k != -1, "Parameter k not set");
            check_true(k+1 <= KMER_MAX_LENGTH, "Maximum allowed k is " + std::to_string(KMER_MAX_LENGTH - 1) + ". To increase the limit, recompile by first running cmake with the option `-DMAX_KMER_LENGTH=n`, where n is a number up to 255, and then running `make` again."); // 255 is max because of KMC
        }

        check_writable(index_dbg_file);
        check_writable(index_color_file);

        if(colorfile != ""){
            check_true(!no_colors, "Must not give both --no-colors and --colorfile");
            check_readable(colorfile);
        }

        check_true(temp_dir != "", "Temp directory not set");
        check_dir_exists(temp_dir);

        check_true(memory_megas > 0, "Memory budget must be positive");
        check_true(colorset_sampling_distance >= 1, "Colorset sampling distance must be positive");

    }

    string to_string(){
        stringstream ss;
        ss << "Input file = " << inputfile << "\n";
        ss << "Input format = " << input_format << "\n";
        if(colorfile != "") ss << "Color name file = " << colorfile << "\n";
        ss << "Index de Bruijn graph output file = " << index_dbg_file << "\n";
        ss << "Index coloring output file = " << index_color_file << "\n";
        ss << "Temporary directory = " << temp_dir << "\n";
        ss << "k = " << k << "\n";
        ss << "Number of threads = " << n_threads << "\n";
        ss << "Memory megabytes = " << memory_megas << "\n";
        ss << "User-specified colors = " << (colorfile == "" ? "false" : "true") << "\n";
        ss << "Load DBG = " << (load_dbg ? "true" : "false") << "\n";
        ss << "Handing of non-ACGT characters = " << (del_non_ACGT ? "delete" : "randomize"); // Last has no endline
        return ss.str();
    }
};

int build_index_main(int argc, char** argv){

    // Legacy support: transform old option format --k to -k
    string legacy_support_fix = "-k";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--auto-colors"){
            cerr << "Flag --auto-colors is now redundant as it is the default behavior starting from Themisto 2.0.0. There are also other changes to the CLI in Themisto 2.0.0. Please read the release notes and update your scripts accordingly." << endl;
            return 1;
        }
        if(string(argv[i]) == "--k") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "Builds an index consisting of compact de Bruijn graph using the Wheeler graph data structure and color information. The input is a set of reference sequences in a single file in fasta or fastq format, and a colorfile, which is a plain text file containing the colors (integers) of the reference sequences in the same order as they appear in the reference sequence file, one line per sequence. If there are characters outside of the DNA alphabet ACGT in the input sequences, those are replaced with random characters from the DNA alphabet.");

    options.add_options()
        ("k,node-length", "The k of the k-mers.", cxxopts::value<LL>())
        ("i,input-file", "The input sequences in FASTA or FASTQ format. The format is inferred from the file extension. Recognized file extensions for fasta are: .fasta, .fna, .ffn, .faa and .frn . Recognized extensions for fastq are: .fastq and .fq . If the file ends with .gz, it is uncompressed into a temporary directory and the temporary file is deleted after use.", cxxopts::value<string>())
        ("c,color-file", "One color per sequence in the fasta file, one color per line. If not given, the sequences are given colors 0,1,2... in the order they appear in the input file.", cxxopts::value<string>()->default_value(""))
        ("o,index-prefix", "The de Bruijn graph will be written to [prefix].themisto.dbg and the color structure to [prefix].themisto.colors.", cxxopts::value<string>())
        ("temp-dir", "Directory for temporary files. This directory should have fast I/O operations and should have as much space as possible.", cxxopts::value<string>())
        ("m,mem-megas", "Number of megabytes allowed for external memory algorithms. Default: 1000", cxxopts::value<LL>()->default_value("1000"))
        ("t,n-threads", "Number of parallel exectuion threads. Default: 1", cxxopts::value<LL>()->default_value("1"))
        ("randomize-non-ACGT", "Replace non-ACGT letters with random nucleotides. If this option is not given, (k+1)-mers containing a non-ACGT character are deleted instead.", cxxopts::value<bool>()->default_value("false"))
        ("d,colorset-pointer-tradeoff", "This option controls a time-space tradeoff for storing and querying color sets. If given a value d, we store color set pointers only for every d nodes on every unitig. The higher the value of d, the smaller then index, but the slower the queries. The savings might be significant if the number of distinct color sets is small and the graph is large and has long unitigs.", cxxopts::value<LL>()->default_value("1"))
        ("no-colors", "Build only the de Bruijn graph without colors.", cxxopts::value<bool>()->default_value("false"))
        ("load-dbg", "If given, loads a precomputed de Bruijn graph from the index directory. If this is given, the parameter -k must not be given because the order k is defined by the precomputed de Bruijn graph.", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        cerr << "Usage examples:" << endl;
        cerr << "Build the de Bruijn graph and colors:" << endl;
        cerr << "  " << argv[0] << " -k 31 --mem-megas 10000 --input-file references.fna --color-file colors.txt --index-prefix my_index --temp-dir temp" << endl;
        cerr << "Build only the de Bruijn graph" << endl;
        cerr << "  " << argv[0] << " -k 31 --mem-megas 10000 --input-file references.fna --index-prefix my_index --temp-dir temp --no-colors" << endl;
        cerr << "Load a previously built de Bruijn graph from the index directory and compute the colors:" << endl;
        cerr << "  " << argv[0] << " --mem-megas 10000 --input-file references.fna --color-file colors.txt --index-prefix my_index --temp-dir temp --load-dbg" << endl;
        exit(1);
    }    

    Build_Config C;
    Themisto themisto;
    C.k = opts["k"].as<LL>();
    C.inputfile = opts["input-file"].as<string>();
    C.input_format = figure_out_file_format(C.inputfile);
    C.n_threads = opts["n-threads"].as<LL>();
    C.colorfile = opts["color-file"].as<string>();
    C.index_dbg_file = opts["index-prefix"].as<string>() + ".themisto.dbg";
    C.index_color_file = opts["index-prefix"].as<string>() + ".themisto.colors";
    C.temp_dir = opts["temp-dir"].as<string>();
    C.load_dbg = opts["load-dbg"].as<bool>();
    C.memory_megas = opts["mem-megas"].as<LL>();
    C.no_colors = opts["no-colors"].as<bool>();
    C.colorset_sampling_distance = opts["colorset-pointer-tradeoff"].as<LL>();
    C.del_non_ACGT = !(opts["randomize-non-ACGT"].as<bool>());

    create_directory_if_does_not_exist(C.temp_dir);

    C.check_valid();
    get_temp_file_manager().set_dir(C.temp_dir);

    cerr << C.to_string() << endl;
    write_log("Starting");

    if(C.input_format == "gzip"){
        write_log("Decompressing the input file");
        string new_name = get_temp_file_manager().create_filename("input");
        check_true(gz_decompress(C.inputfile, new_name) == Z_OK, "Problem with zlib decompression");
        C.input_format = figure_out_file_format(C.inputfile.substr(0,C.inputfile.size() - 3));
        C.inputfile = new_name;
    }

    if(!C.no_colors && C.colorfile == ""){
        // Automatic colors
        write_log("Assigning colors");
        C.colorfile = generate_default_colorfile(C.inputfile, C.input_format);
    }

    // Deal with non-ACGT characters
    if(C.del_non_ACGT){
        std::tie(C.inputfile, C.colorfile) = split_all_seqs_at_non_ACGT(C.inputfile, C.input_format, C.colorfile); // Turns the file into fasta format also
        C.input_format = "fasta"; // split_all_seqs_at_non_ACGT returns a fasta file
    } else {
        C.inputfile = fix_alphabet(C.inputfile, C.input_format == "fasta" ? FASTA_MODE : FASTQ_MODE); // Turns the file into fasta format also
        C.input_format = "fasta"; // fix_alphabet returns a fasta file
    }
    
    if(C.load_dbg){
        write_log("Loading de Bruijn Graph");
        themisto.load_boss(C.index_dbg_file);
    } else{
        write_log("Building de Bruijn Graph");
        themisto.construct_boss(C.inputfile, C.k, C.memory_megas * 1e6, C.n_threads, false);
        themisto.save_boss(C.index_dbg_file);
        write_log("Building de Bruijn Graph finished (" + std::to_string(themisto.boss.number_of_nodes()) + " nodes)");
    }

    if(!C.no_colors){
        write_log("Building colors");
        themisto.construct_colors(C.inputfile, C.colorfile, C.memory_megas * 1e6, C.n_threads, C.colorset_sampling_distance);
        themisto.save_colors(C.index_color_file);
    } else{
        std::filesystem::remove(C.index_color_file); // There is an empty file so let's remove it
    }

    write_log("Finished");

    return 0;
}