//#include "Themisto.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "sbwt/globals.hh"
#include "sbwt/variants.hh"
#include "new_coloring.hh"
#include "globals.hh"
#include "sbwt/SeqIO.hh"
#include "sbwt/cxxopts.hpp"

using namespace std;
typedef long long LL;

struct Build_Config{
    LL k = 0;
    LL n_threads = 1;
    string inputfile;
    string colorfile;
    string index_dbg_file;
    string index_color_file;
    string temp_dir;
    string coloring_structure_type;
    sbwt::SeqIO::FileFormat input_format;
    bool load_dbg = false;
    LL memory_megas = 1000;
    bool no_colors = false;
    bool del_non_ACGT = false;
    LL colorset_sampling_distance = 1;
    bool verbose = false;
    bool silent = false;
    
    void check_valid(){
        sbwt::check_true(inputfile != "", "Input file not set");

        sbwt::check_readable(inputfile);

        if(!load_dbg){
            sbwt::check_true(k != 0, "Parameter k not set");
            sbwt::check_true(k+1 <= KMER_MAX_LENGTH, "Maximum allowed k is " + std::to_string(KMER_MAX_LENGTH - 1) + ". To increase the limit, recompile by first running cmake with the option `-DMAX_KMER_LENGTH=n`, where n is a number up to 255, and then running `make` again."); // 255 is max because of KMC
        } else{
            if(k != 0){
                sbwt::write_log("Warning: value of parameter k is ignored because the DBG is not built, but loaded from disk instead", sbwt::LogLevel::MAJOR);
            }
        }

        sbwt::check_writable(index_dbg_file);
        sbwt::check_writable(index_color_file);

        if(colorfile != ""){
            sbwt::check_true(!no_colors, "Must not give both --no-colors and --colorfile");
            sbwt::check_readable(colorfile);
        }

        if(coloring_structure_type != "sdsl-fixed" && coloring_structure_type != "sdsl-delta" && coloring_structure_type != "roaring"){
            throw std::runtime_error("Unknown coloring structure type: " + coloring_structure_type);
        }

        sbwt::check_true(temp_dir != "", "Temp directory not set");
        check_dir_exists(temp_dir);

        sbwt::check_true(memory_megas > 0, "Memory budget must be positive");
        sbwt::check_true(colorset_sampling_distance >= 1, "Colorset sampling distance must be positive");

    }

    string to_string(){
        stringstream ss;
        ss << "Input file = " << inputfile << "\n";
        ss << "Input format = " << input_format.extension << "\n";
        if(colorfile != "") ss << "Color name file = " << colorfile << "\n";
        ss << "Index de Bruijn graph output file = " << index_dbg_file << "\n";
        ss << "Index coloring output file = " << index_color_file << "\n";
        ss << "Temporary directory = " << temp_dir << "\n";
        ss << "k = " << k << "\n";
        ss << "Number of threads = " << n_threads << "\n";
        ss << "Memory megabytes = " << memory_megas << "\n";
        ss << "User-specified colors = " << (colorfile == "" ? "false" : "true") << "\n";
        ss << "Load DBG = " << (load_dbg ? "true" : "false") << "\n";
        ss << "Handling of non-ACGT characters = " << (del_non_ACGT ? "delete" : "randomize") << "\n";
        ss << "Coloring structure type: " << coloring_structure_type << "\n"; 

        string verbose_level = "normal";
        if(verbose) verbose_level = "verbose";
        if(silent) verbose_level = "silent";

        ss << "Verbosity = " << verbose_level; // Last has no endline

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
        ("k,node-length", "The k of the k-mers.", cxxopts::value<LL>()->default_value("0"))
        ("i,input-file", "The input sequences in FASTA or FASTQ format. The format is inferred from the file extension. Recognized file extensions for fasta are: .fasta, .fna, .ffn, .faa and .frn . Recognized extensions for fastq are: .fastq and .fq . If the file ends with .gz, it is uncompressed into a temporary directory and the temporary file is deleted after use.", cxxopts::value<string>())
        ("c,color-file", "One color per sequence in the fasta file, one color per line. If not given, the sequences are given colors 0,1,2... in the order they appear in the input file.", cxxopts::value<string>()->default_value(""))
        ("o,index-prefix", "The de Bruijn graph will be written to [prefix].tdbg and the color structure to [prefix].tcolors.", cxxopts::value<string>())
        ("temp-dir", "Directory for temporary files. This directory should have fast I/O operations and should have as much space as possible.", cxxopts::value<string>())
        ("m,mem-megas", "Number of megabytes allowed for external memory algorithms (must be at least 2048).", cxxopts::value<LL>()->default_value("2048"))
        ("t,n-threads", "Number of parallel exectuion threads. Default: 1", cxxopts::value<LL>()->default_value("1"))
        ("randomize-non-ACGT", "Replace non-ACGT letters with random nucleotides. If this option is not given, (k+1)-mers containing a non-ACGT character are deleted instead.", cxxopts::value<bool>()->default_value("false"))
        ("d,colorset-pointer-tradeoff", "This option controls a time-space tradeoff for storing and querying color sets. If given a value d, we store color set pointers only for every d nodes on every unitig. The higher the value of d, the smaller then index, but the slower the queries. The savings might be significant if the number of distinct color sets is small and the graph is large and has long unitigs.", cxxopts::value<LL>()->default_value("1"))
        ("no-colors", "Build only the de Bruijn graph without colors.", cxxopts::value<bool>()->default_value("false"))
        ("load-dbg", "If given, loads a precomputed de Bruijn graph from the index prefix. If this is given, the value of parameter -k is ignored because the order k is defined by the precomputed de Bruijn graph.", cxxopts::value<bool>()->default_value("false"))
        ("s,coloring-structure-type", "Type of coloring structure to build (\"sdsl-fixed\", \"sdsl-delta\" or \"roaring\" ).", cxxopts::value<string>()->default_value("sdsl-delta"))
        ("v,verbose", "More verbose progress reporting into stderr.", cxxopts::value<bool>()->default_value("false"))
        ("silent", "Print as little as possible to stderr (only errors).", cxxopts::value<bool>()->default_value("false"))
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
    C.k = opts["k"].as<LL>();
    C.inputfile = opts["input-file"].as<string>();
    C.input_format = sbwt::SeqIO::figure_out_file_format(C.inputfile);
    C.n_threads = opts["n-threads"].as<LL>();
    C.colorfile = opts["color-file"].as<string>();
    C.index_dbg_file = opts["index-prefix"].as<string>() + ".tdbg";
    C.index_color_file = opts["index-prefix"].as<string>() + ".tcolors";
    C.temp_dir = opts["temp-dir"].as<string>();
    C.load_dbg = opts["load-dbg"].as<bool>();
    C.memory_megas = opts["mem-megas"].as<LL>();
    C.no_colors = opts["no-colors"].as<bool>();
    C.colorset_sampling_distance = opts["colorset-pointer-tradeoff"].as<LL>();
    C.del_non_ACGT = !(opts["randomize-non-ACGT"].as<bool>());
    C.verbose = opts["verbose"].as<bool>();
    C.silent = opts["silent"].as<bool>();
    C.coloring_structure_type = opts["coloring-structure-type"].as<string>();

    if(C.verbose && C.silent) throw runtime_error("Can not give both --verbose and --silent");
    if(C.verbose) set_log_level(sbwt::LogLevel::MINOR);
    if(C.silent) set_log_level(sbwt::LogLevel::OFF);

    create_directory_if_does_not_exist(C.temp_dir);

    C.check_valid();
    sbwt::get_temp_file_manager().set_dir(C.temp_dir);

    write_log("Build configuration:\n" + C.to_string(), sbwt::LogLevel::MAJOR);
    write_log("Starting", sbwt::LogLevel::MAJOR);

    // Deal with non-ACGT characters
    if(C.del_non_ACGT){
        // KMC takes care of this
    } else {
        write_log("Replacing non-ACGT characters with random nucleotides", LogLevel::MAJOR);
        C.inputfile = fix_alphabet(C.inputfile); // Turns the file into fasta format also
    }

    if(!C.no_colors && C.colorfile == ""){
        // Automatic colors
        sbwt::write_log("Assigning colors", sbwt::LogLevel::MAJOR);
        C.colorfile = generate_default_colorfile(C.inputfile);
    }

    std::unique_ptr<sbwt::plain_matrix_sbwt_t> dbg_ptr;

    if(C.load_dbg){
        sbwt::write_log("Loading de Bruijn Graph", sbwt::LogLevel::MAJOR);
        dbg_ptr = std::make_unique<sbwt::plain_matrix_sbwt_t>();
        dbg_ptr->load(C.index_dbg_file);
    } else{
        sbwt::write_log("Building de Bruijn Graph", sbwt::LogLevel::MAJOR);
        sbwt::plain_matrix_sbwt_t::BuildConfig sbwt_config;
        sbwt_config.add_reverse_complements = false;
        sbwt_config.build_streaming_support = true;
        sbwt_config.input_files = {C.inputfile};
        sbwt_config.k = C.k;
        sbwt_config.max_abundance = 1e9;
        sbwt_config.min_abundance = 1;
        sbwt_config.n_threads = C.n_threads;
        sbwt_config.ram_gigas = min(2LL, C.memory_megas / (1 << 10)); // KMC requires at least 2 GB
        sbwt_config.temp_dir = C.temp_dir;
        dbg_ptr = std::make_unique<sbwt::plain_matrix_sbwt_t>(sbwt_config);
        dbg_ptr->serialize(C.index_dbg_file);
        sbwt::write_log("Building de Bruijn Graph finished (" + std::to_string(dbg_ptr->number_of_kmers()) + " k-mers)", sbwt::LogLevel::MAJOR);
    }

    if(!C.no_colors){
        sbwt::write_log("Building colors", sbwt::LogLevel::MAJOR);

        vector<int64_t> color_assignment = read_colorfile(C.colorfile);
        if(C.coloring_structure_type == "sdsl-fixed"){
            Coloring<Fixed_Width_Int_Color_Set> coloring;
            coloring.add_colors(*dbg_ptr, C.inputfile, color_assignment, C.memory_megas * (1 << 20), C.n_threads, C.colorset_sampling_distance);
            sbwt::throwing_ofstream out(C.index_color_file, ios::binary);
            coloring.serialize(out.stream);
        } if(C.coloring_structure_type == "sdsl-delta"){
            Coloring<Bitmap_Or_Deltas_ColorSet> coloring;
            coloring.add_colors(*dbg_ptr, C.inputfile, color_assignment, C.memory_megas * (1 << 20), C.n_threads, C.colorset_sampling_distance);
            sbwt::throwing_ofstream out(C.index_color_file, ios::binary);
            coloring.serialize(out.stream);
        } else if(C.coloring_structure_type == "roaring"){
            Coloring<Roaring_Color_Set> coloring;
            coloring.add_colors(*dbg_ptr, C.inputfile, color_assignment, C.memory_megas * (1 << 20), C.n_threads, C.colorset_sampling_distance);
            sbwt::throwing_ofstream out(C.index_color_file, ios::binary);
            coloring.serialize(out.stream);
        }
    } else{
        std::filesystem::remove(C.index_color_file); // There is an empty file so let's remove it
    }

    sbwt::write_log("Finished", sbwt::LogLevel::MAJOR);

    return 0;
}