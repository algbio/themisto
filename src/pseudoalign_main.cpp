#include "Themisto.hh"
#include "input_reading.hh"
#include "zpipe.hh"
#include <string>
#include <cstring>
#include "version.h"
#include "cxxopts.hpp"

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
        for(string query_file : query_files){
            if(query_file != ""){
                check_readable(query_file);
            }
        }

        for(string outfile : outfiles){
            check_true(outfile != "", "Outfile not set");
            check_writable(outfile);
        }

        check_true(query_files.size() == outfiles.size(), "Number of query files and outfiles do not match");
        
        check_true(index_dir != "", "Index directory not set");
        check_dir_exists(index_dir);

        check_true(temp_dir != "", "Temp directory not set");
        check_dir_exists(temp_dir);
    }
};

char get_rc(char c){
    switch(c){
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return c;
    }
}

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
        ("o,out-file", "Output filename.", cxxopts::value<string>()->default_value(""))
        ("out-file-list", "A file containing a list of output filenames, one per line.", cxxopts::value<string>()->default_value(""))
        ("i,index-dir", "Directory where the index will be built. Always required, directory must exist before running.", cxxopts::value<string>())
        ("temp-dir", "Directory for temporary files.", cxxopts::value<string>())
        ("rc", "Whether to to consider the reverse complement k-mers in the pseudoalignemt.", cxxopts::value<bool>()->default_value("false")) 
        ("t, n-threads", "Number of parallel exectuion threads. Default: 1", cxxopts::value<LL>()->default_value("1"))
        ("gzip-output", "Compress the output files with gzip.", cxxopts::value<bool>()->default_value("false"))
        ("sort-output", "Sort the lines of the out files by sequence rank in the input files.", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    LL old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        cerr << "Usage examples:" << endl;
        cerr << "Pseudoalign reads.fna against an index:" << endl;
        cerr << "  pseudoalign --query-file reads.fna --index-dir index --temp-dir temp --out-file out.txt" << endl;
        cerr << "Pseudoalign a list of fasta files in input_list.txt into output filenames in output_list.txt:" << endl;
        cerr << "  pseudoalign --query-file-list input_list.txt --index-dir index --temp-dir temp --out-file-list output_list.txt" << endl;
        cerr << "Pseudoalign reads.fna against an index using also reverse complements:" << endl;
        cerr << "  pseudoalign --rc --query-file reads.fna --index-dir index --temp-dir temp --outfile out.txt" << endl;
        exit(1);
    }

    Config C;
    if(opts.count("query-file") && opts["query-file"].as<string>() != "") C.query_files.push_back(opts["query-file"].as<string>());
    if(opts.count("query-file-list") && opts["query-file-list"].as<string>() != "") 
        for(string line : read_lines(opts["query-file-list"].as<string>()))
            C.query_files.push_back(line);
    if(opts.count("out-file") && opts["out-file"].as<string>() != "") C.outfiles.push_back(opts["out-file"].as<string>());
    if(opts.count("out-file-list") && opts["out-file-list"].as<string>() != "") 
        for(string line : read_lines(opts["out-file-list"].as<string>()))
            C.outfiles.push_back(line);
    C.index_dir = opts["index-dir"].as<string>();
    C.temp_dir = opts["temp-dir"].as<string>();
    C.reverse_complements = opts["rc"].as<bool>();
    C.n_threads = opts["n-threads"].as<LL>();
    C.gzipped_output = opts["gzip-output"].as<bool>();
    C.sort_output = opts["sort-output"].as<bool>();

    create_directory_if_does_not_exist(C.temp_dir);

    if(C.gzipped_output){
        for(string& filename : C.outfiles) filename += ".gz";
    }

    C.check_valid();

    write_log("Starting");

    temp_file_manager.set_dir(C.temp_dir);

    write_log("Loading the index");    
    Themisto themisto;
    themisto.load_boss(C.index_dir + "/boss-");
    themisto.load_colors(C.index_dir + "/coloring-");

    for(LL i = 0; i < C.query_files.size(); i++){
        write_log("Aligning " + C.query_files[i] + " (writing output to " + C.outfiles[i] + ")");

        string inputfile = C.query_files[i];
        string file_format = figure_out_file_format(inputfile);
        if(file_format == "gzip"){
            string new_name = temp_file_manager.get_temp_file_name("input");
            check_true(gz_decompress(inputfile, new_name) == Z_OK, "Problem with zlib decompression");
            file_format = figure_out_file_format(inputfile.substr(0,inputfile.size() - 3));
            inputfile = new_name;
        }

        Sequence_Reader sr(inputfile, file_format == "fasta" ? FASTA_MODE : FASTQ_MODE);
        sr.set_upper_case(true);
        themisto.pseudoalign_parallel(C.n_threads, sr, C.outfiles[i], C.reverse_complements, 1000000, C.gzipped_output, C.sort_output); // Buffer size 1 MB
        temp_file_manager.clean_up();
    }

    write_log("Finished");

    return 0;
}
