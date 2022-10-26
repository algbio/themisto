#include <string>
#include <cstring>
#include "zpipe.hh"
#include "version.h"
#include "Coloring.hh"
#include "globals.hh"
#include "pseudoalign.hh"
#include "sbwt/globals.hh"
#include "sbwt/throwing_streams.hh"
#include "sbwt/variants.hh"
#include "sbwt/SeqIO.hh"
#include "sbwt/cxxopts.hpp"

using namespace std;
using namespace sbwt;

typedef long long LL;

struct Pseudoalign_Config{
    vector<string> query_files;
    vector<string> outfiles;
    string index_dbg_file;
    string index_color_file;
    string temp_dir;

    bool gzipped_output = false;
    bool reverse_complements = false;
    bool sort_output = false;
    LL n_threads = 1;
    double buffer_size_megas = 8;
    bool verbose = false;
    bool silent = false;

    void check_valid(){
        for(string query_file : query_files){
            if(query_file != ""){
                check_readable(query_file);
            }
        }

        check_readable(index_dbg_file);
        check_readable(index_color_file);

        for(string outfile : outfiles){
            check_true(outfile != "", "Outfile not set");
            check_writable(outfile);
        }

    if (outfiles.size() == 0 && query_files.size() > 1) {
        check_true(query_files.size() == 1, "Can't print results when aligning multiple files; supply " + to_string(query_files.size()) + " outfiles with --out-file-list.");
    }
    if (outfiles.size() > 0) {
        check_true(query_files.size() == outfiles.size(), "Number of query files and outfiles do not match");
    }

    if (sort_output) {
        check_true(outfiles.size() > 0, "Can't sort output when printing results");
    }
    check_true(temp_dir != "", "Temp directory not set");
        check_dir_exists(temp_dir);
    }
};

vector<string> read_lines(string filename){
    check_readable(filename);
    vector<string> lines;
    sbwt::throwing_ifstream in(filename);
    string line;
    while(in.getline(line)){
        lines.push_back(line);
    }
    return lines;
}

// If outputfile is an empty string, prints to stdout
template<typename coloring_t> 
void call_pseudoalign(plain_matrix_sbwt_t& SBWT, const coloring_t& coloring, Pseudoalign_Config& C, string inputfile, string outputfile){
    if(SeqIO::figure_out_file_format(inputfile).gzipped){
        SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> reader(inputfile);
        pseudoalign(SBWT, coloring, C.n_threads, reader, outputfile, C.reverse_complements, C.buffer_size_megas * (1 << 20), C.gzipped_output, C.sort_output); // Buffer size 8 MB
    } else{
        SeqIO::Reader<Buffered_ifstream<std::ifstream>> reader(inputfile);
        pseudoalign(SBWT, coloring, C.n_threads, reader, outputfile, C.reverse_complements, C.buffer_size_megas * (1 << 20), C.gzipped_output, C.sort_output); // Buffer size 8 MB
    }
}

int pseudoalign_main(int argc, char** argv){

    // Legacy support: transform old option format --outfile --out-file
    string legacy_support_fix = "--out-file";
    for(LL i = 1; i < argc; i++){
        if(string(argv[i]) == "--outfile") argv[i] = &(legacy_support_fix[0]);
    }

    cxxopts::Options options(argv[0], "This program aligns query sequences against an index that has been built previously. The output is one line per input read. Each line consists of a space-separated list of integers. The first integer specifies the rank of the read in the input file, and the rest of the integers are the identifiers of the colors of the sequences that the read pseudoaligns with. If the program is ran with more than one thread, the output lines are not necessarily in the same order as the reads in the input file. This can be fixed with the option --sort-output.\n\nIf the coloring data structure was built with the --color-file option, then the integer identifiers of the colors can be mapped back to the provided color names by parsing the file coloring-mapping-id_to_name in the index directory. This file contains as many lines as there are distinct colors, and each line contains two space-separated strings: the first is the integer identifier of a color, and the second is the corresponding color name. In case the --auto-colors option was used, the integer identifiers are always numbers [0..n-1], where n is the total number of reference sequences, and the identifiers are assigned in the same order as the reference sequences were given to build_index.\n\n The query can be given as one file, or as a file with a list of files. In the former case, we must specify one output file with the options --out-file, and in the latter case, we must give a file that lists one output filename per line using the option --out-file-list.\n\nThe query file(s) should be in fasta of fastq format. The format is inferred from the file extension. Recognized file extensions for fasta are: .fasta, .fna, .ffn, .faa and .frn . Recognized extensions for fastq are: .fastq and .fq");

    options.add_options()
        ("q, query-file", "Input file of the query sequences", cxxopts::value<string>()->default_value(""))
        ("query-file-list", "A list of query filenames, one line per filename", cxxopts::value<string>()->default_value(""))
        ("o,out-file", "Output filename. Print results if no output filename is given.", cxxopts::value<string>()->default_value(""))
        ("out-file-list", "A file containing a list of output filenames, one per line.", cxxopts::value<string>()->default_value(""))
        ("i,index-prefix", "The index prefix that was given to the build command.", cxxopts::value<string>())
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
        ("rc", "Whether to to consider the reverse complement k-mers in the pseudoalignment.", cxxopts::value<bool>()->default_value("false"))
        ("t, n-threads", "Number of parallel exectuion threads. Default: 1", cxxopts::value<LL>()->default_value("1"))
        ("gzip-output", "Compress the output files with gzip.", cxxopts::value<bool>()->default_value("false"))
        ("sort-output", "Sort the lines of the out files by sequence rank in the input files.", cxxopts::value<bool>()->default_value("false"))
        ("buffer-size-megas", "Size of the input buffer in megabytes in each thread. If this is larger than the number of nucleotides in the input divided by the number of threads, then some threads will be idle. So if your input files are really small and you have a lot of threads, consider using a small buffer.", cxxopts::value<double>()->default_value("8.0"))
        ("v,verbose", "More verbose progress reporting into stderr.", cxxopts::value<bool>()->default_value("false"))
        ("silent", "Print as little as possible to stderr (only errors).", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        cerr << "Usage examples:" << endl;
        cerr << "Pseudoalign reads.fna against an index:" << endl;
        cerr << "  " << argv[0] << " --query-file reads.fna --index-prefix my_index --temp-dir temp --out-file out.txt" << endl;
        cerr << "Pseudoalign a list of fasta files in input_list.txt into output filenames in output_list.txt:" << endl;
        cerr << "  " << argv[0] << " --query-file-list input_list.txt --index-prefix my_index --temp-dir temp --out-file-list output_list.txt" << endl;
        cerr << "Pseudoalign reads.fna against an index using also reverse complements:" << endl;
        cerr << "  " << argv[0] << " --rc --query-file reads.fna --index-prefix my_index --temp-dir temp --outfile out.txt" << endl;
        exit(1);
    }

    Pseudoalign_Config C;
    if(opts.count("query-file") && opts["query-file"].as<string>() != "") C.query_files.push_back(opts["query-file"].as<string>());
    if(opts.count("query-file-list") && opts["query-file-list"].as<string>() != "")
        for(string line : read_lines(opts["query-file-list"].as<string>()))
            C.query_files.push_back(line);
    if(opts.count("out-file") && opts["out-file"].as<string>() != "") C.outfiles.push_back(opts["out-file"].as<string>());
    if(opts.count("out-file-list") && opts["out-file-list"].as<string>() != "")
        for(string line : read_lines(opts["out-file-list"].as<string>()))
            C.outfiles.push_back(line);
    C.index_dbg_file = opts["index-prefix"].as<string>() + ".tdbg";
    C.index_color_file = opts["index-prefix"].as<string>() + ".tcolors";
    C.temp_dir = opts["temp-dir"].as<string>();
    C.reverse_complements = opts["rc"].as<bool>();
    C.n_threads = opts["n-threads"].as<LL>();
    C.gzipped_output = opts["gzip-output"].as<bool>();
    C.sort_output = opts["sort-output"].as<bool>();
    C.verbose = opts["verbose"].as<bool>();
    C.silent = opts["silent"].as<bool>();
    C.buffer_size_megas = opts["buffer-size-megas"].as<double>();

    if(C.verbose && C.silent) throw runtime_error("Can not give both --verbose and --silent");
    if(C.verbose) set_log_level(LogLevel::MINOR);
    if(C.silent) set_log_level(LogLevel::OFF);

    create_directory_if_does_not_exist(C.temp_dir);

    if(C.gzipped_output){
        for(string& filename : C.outfiles) filename += ".gz";
    }

    C.check_valid();

    write_log("Starting", LogLevel::MAJOR);

    get_temp_file_manager().set_dir(C.temp_dir);

    write_log("Loading the index", LogLevel::MAJOR);
    plain_matrix_sbwt_t SBWT;
    SBWT.load(C.index_dbg_file);

    // Load whichever coloring data structure type is stored on disk
    std::variant<Coloring<SDSL_Variant_Color_Set>,
                 Coloring<Roaring_Color_Set>,
                 Coloring<Bit_Magic_Color_Set>> coloring;
    load_coloring(C.index_color_file, SBWT, coloring);

    if(std::holds_alternative<Coloring<SDSL_Variant_Color_Set>>(coloring))
        write_log("sdsl coloring structure loaded", LogLevel::MAJOR);
    if(std::holds_alternative<Coloring<Roaring_Color_Set>>(coloring))
        write_log("roaring coloring structure loaded", LogLevel::MAJOR);
    if(std::holds_alternative<Coloring<Bit_Magic_Color_Set>>(coloring))
        write_log("bitmagic coloring structure loaded", LogLevel::MAJOR);

    for(LL i = 0; i < C.query_files.size(); i++){
        if (C.outfiles.size() > 0) {
            write_log("Aligning " + C.query_files[i] + " (writing output to " + C.outfiles[i] + ")", LogLevel::MAJOR);
        } else {
            write_log("Aligning " + C.query_files[i] + " (printing output)", LogLevel::MAJOR);
        }

        if(std::holds_alternative<Coloring<SDSL_Variant_Color_Set>>(coloring))
            call_pseudoalign(SBWT, get<Coloring<SDSL_Variant_Color_Set>>(coloring), C, C.query_files[i], (C.outfiles.size() > 0 ? C.outfiles[i] : ""));
        if(std::holds_alternative<Coloring<Roaring_Color_Set>>(coloring))
            call_pseudoalign(SBWT, get<Coloring<Roaring_Color_Set>>(coloring), C, C.query_files[i], (C.outfiles.size() > 0 ? C.outfiles[i] : ""));
        if(std::holds_alternative<Coloring<Bit_Magic_Color_Set>>(coloring))
            call_pseudoalign(SBWT, get<Coloring<Bit_Magic_Color_Set>>(coloring), C, C.query_files[i], (C.outfiles.size() > 0 ? C.outfiles[i] : ""));
    }

    write_log("Finished", LogLevel::MAJOR);

    return 0;
}
