#include <string>
#include <cstring>
#include "zpipe.hh"
#include "version.h"
#include "coloring/Coloring.hh"
#include "globals.hh"
#include "pseudoalign.hh"
#include "sbwt/globals.hh"
#include "sbwt/throwing_streams.hh"
#include "sbwt/variants.hh"
#include "SeqIO/SeqIO.hh"
#include "sbwt/cxxopts.hpp"

using namespace std;
using namespace sbwt;

struct Pseudoalign_Config{
    vector<string> query_files;
    vector<string> outfiles;
    string index_dbg_file;
    string index_color_file;
    string temp_dir;
    string aux_info_file;

    bool gzipped_output = false;
    bool reverse_complements = false;
    bool sort_output_lines = false;
    bool sort_hits = false;
    int64_t n_threads = 1;
    double buffer_size_megas = 8;
    bool verbose = false;
    bool silent = false;
    double threshold = -1;
    bool ignore_unknown = false;
    double relevant_kmers_fraction = 0;

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

    if (sort_output_lines) {
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
    if(seq_io::figure_out_file_format(inputfile).gzipped){
        seq_io::Reader<seq_io::Buffered_ifstream<seq_io::zstr::ifstream>> reader(inputfile);
        pseudoalign(SBWT, coloring, C.n_threads, reader, outputfile, C.aux_info_file, C.reverse_complements, C.buffer_size_megas * (1 << 20), C.gzipped_output, C.sort_output_lines, C.threshold, C.ignore_unknown, C.relevant_kmers_fraction, C.sort_hits); // Buffer size 8 MB
    } else{
        seq_io::Reader<seq_io::Buffered_ifstream<std::ifstream>> reader(inputfile);
        pseudoalign(SBWT, coloring, C.n_threads, reader, outputfile, C.aux_info_file, C.reverse_complements, C.buffer_size_megas * (1 << 20), C.gzipped_output, C.sort_output_lines, C.threshold, C.ignore_unknown, C.relevant_kmers_fraction, C.sort_hits); // Buffer size 8 MB
    }
}

int pseudoalign_main(int argc_given, char** argv_given){

    // Legacy support: transform old options
    char** argv = (char**)malloc(sizeof(char*) * argc_given); // Freed and the and of the function
    argv[0] = argv_given[0];
    int64_t argc = 1;
    char legacy_support_fix[] = "--out-file";
    char legacy_support_fix2[] = "--sort-output-lines"; // --sort-output is now this
    for(int64_t i = 1; i < argc_given; i++){
        if(string(argv_given[i]) == "--outfile") argv[argc++] = legacy_support_fix;
        else if(string(argv_given[i]) == "--sort-output"){
            cerr << "Note: --sort-output is now called " + string(legacy_support_fix2) << ", but the old option name still works." <<  endl;
            argv[argc++] = legacy_support_fix2;
        }
        else if(string(argv_given[i]) == "--ignore-unknown-kmers"){
            // This is now the default. Remove (ignore) the flag.
        } else argv[argc++] = argv_given[i];
    }

    cxxopts::Options options(argv[0], "This program aligns query sequences against an index that has been built previously. The output is one line per input read. Each line consists of a space-separated list of integers. The first integer specifies the rank of the read in the input file, and the rest of the integers are the identifiers of the colors of the sequences that the read pseudoaligns with. If the program is ran with more than one thread, the output lines are not necessarily in the same order as the reads in the input file. This can be fixed with the option --sort-output, but this will slow down the program.\n\n The query can be given as one file, or as a file with a list of files. In the former case, we must specify one output file with the options --out-file, and in the latter case, we must give a file that lists one output filename per line using the option --out-file-list.\n\nThe query file(s) should be in fasta of fastq format. The format is inferred from the file extension. Recognized file extensions for fasta are: .fasta, .fna, .ffn, .faa and .frn . Recognized extensions for fastq are: .fastq and .fq. Gzipped sequence files with the extension .gz are also supported.");

    options.add_options("Basic")
        ("q, query-file", "Input file of the query sequences", cxxopts::value<string>()->default_value(""))
        ("query-file-list", "A list of query filenames, one line per filename", cxxopts::value<string>()->default_value(""))
        ("o,out-file", "Output filename. Print results if no output filename is given.", cxxopts::value<string>()->default_value(""))
        ("out-file-list", "A file containing a list of output filenames, one per line.", cxxopts::value<string>()->default_value(""))
        ("auxiliary-info-file", "Optional: Write to this file auxiliary information for each read. On each line, three space-separated integers: read rank, number of relevant k-mers and number of k-mers.", cxxopts::value<string>())
        ("i,index-prefix", "The index prefix that was given to the build command.", cxxopts::value<string>())
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
        ("gzip-output", "Compress the output files with gzip.", cxxopts::value<bool>()->default_value("false"))
        ("sort-output-lines", "Sort the lines in the output files by sequence rank in the input files. To sort the color ids *within* the lines, use --sort-hits.", cxxopts::value<bool>()->default_value("false"))
        ("sort-hits", "Sort the color ids within each line of the output.", cxxopts::value<bool>()->default_value("false"))
        ("v,verbose", "More verbose progress reporting into stderr.", cxxopts::value<bool>()->default_value("false"))
    ;

    options.add_options("Algorithm")
        ("threshold", "Fraction of k-mer matches required to report a color. If this is equal to 1, the algorithm is implemented with a specialized set intersection method.", cxxopts::value<double>()->default_value("1"))
        ("include-unknown-kmers", "Include all k-mers in the pseudoalignment, even those which do not occur in the index.", cxxopts::value<bool>()->default_value("false"))
        ("relevant-kmers-fraction", "Accept a pseudoalignment only if at least this fraction of k-mers of the read had at least 1 color.", cxxopts::value<double>()->default_value("0.0"))
    ;

    options.add_options("Computational resources")
        ("t, n-threads", "Number of parallel execution threads. Default: 1", cxxopts::value<int64_t>()->default_value("1"))
    ;

    options.add_options("Help")
        ("h,help", "Print usage instructions for commonly used options.")
        ("help-advanced", "Print advanced usage instructions.")
    ;

    options.add_options("Advanced")
        ("rc", "Include reverse complement matches in the pseudoalignment. This option only makes sense if the index was built with --forward-strand-only. Otherwise this option has no effect except to slow down the query.", cxxopts::value<bool>()->default_value("false"))
        ("buffer-size-megas", "Size of the input buffer in megabytes in each thread. If this is larger than the number of nucleotides in the input divided by the number of threads, then some threads will be idle. So if your input files are really small and you have a lot of threads, consider using a small buffer.", cxxopts::value<double>()->default_value("8.0"))
        ("silent", "Print as little as possible to stderr (only errors).", cxxopts::value<bool>()->default_value("false"))
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help") || opts.count("help-advanced")){
        if(old_argc == 1 || opts.count("help"))
            std::cerr << options.help({"Basic","Algorithm","Computational resources","Help"}) << std::endl;
        if(opts.count("help-advanced"))
            std::cerr << options.help({"Basic","Algorithm","Computational resources","Advanced","Help"}) << std::endl;
        cerr << "Usage example:" << endl;
        cerr << argv[0] << " pseudoalign --query-file example_input/queries.fna --index-prefix my_index --temp-dir temp --out-file out.txt --n-threads 4 --threshold 0.7" << endl;
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
    C.n_threads = opts["n-threads"].as<int64_t>();
    C.gzipped_output = opts["gzip-output"].as<bool>();
    C.sort_output_lines = opts["sort-output-lines"].as<bool>();
    C.sort_hits = opts["sort-hits"].as<bool>();
    C.verbose = opts["verbose"].as<bool>();
    C.silent = opts["silent"].as<bool>();
    C.buffer_size_megas = opts["buffer-size-megas"].as<double>();
    C.threshold = opts["threshold"].as<double>();
    C.ignore_unknown = !opts["include-unknown-kmers"].as<bool>();
    C.relevant_kmers_fraction = opts["relevant-kmers-fraction"].as<double>();
    try {
        C.aux_info_file = opts["auxiliary-info-file"].as<string>();
    } catch(cxxopts::option_has_no_value_exception& e){
        // No aux file given -> ok.
    }


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
                 Coloring<Roaring_Color_Set>> coloring;
    load_coloring(C.index_color_file, SBWT, coloring);

    if(std::holds_alternative<Coloring<SDSL_Variant_Color_Set>>(coloring))
        write_log("sdsl coloring structure loaded", LogLevel::MAJOR);
    if(std::holds_alternative<Coloring<Roaring_Color_Set>>(coloring))
        write_log("roaring coloring structure loaded", LogLevel::MAJOR);

    for(int64_t i = 0; i < C.query_files.size(); i++){
        if (C.outfiles.size() > 0) {
            write_log("Aligning " + C.query_files[i] + " (writing output to " + C.outfiles[i] + ")", LogLevel::MAJOR);
        } else {
            write_log("Aligning " + C.query_files[i] + " (printing output)", LogLevel::MAJOR);
        }

        if(std::holds_alternative<Coloring<SDSL_Variant_Color_Set>>(coloring))
            call_pseudoalign(SBWT, get<Coloring<SDSL_Variant_Color_Set>>(coloring), C, C.query_files[i], (C.outfiles.size() > 0 ? C.outfiles[i] : ""));
        if(std::holds_alternative<Coloring<Roaring_Color_Set>>(coloring))
            call_pseudoalign(SBWT, get<Coloring<Roaring_Color_Set>>(coloring), C, C.query_files[i], (C.outfiles.size() > 0 ? C.outfiles[i] : ""));
    }

    write_log("Finished", LogLevel::MAJOR);

    free(argv);
    return 0;
}
