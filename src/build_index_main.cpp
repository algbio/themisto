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
#include "Coloring_Builder.hh"

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
    string from_index;
    sbwt::SeqIO::FileFormat input_format;
    bool load_dbg = false;
    LL memory_megas = 1000;
    bool no_colors = false;
    bool del_non_ACGT = false;
    LL colorset_sampling_distance = 1;
    bool verbose = false;
    bool silent = false;
    bool reverse_complements = false;
    
    void check_valid(){
        sbwt::check_true(inputfile != "", "Input file not set");

        sbwt::check_readable(inputfile);

        if(!load_dbg){
            sbwt::check_true(k != 0, "Parameter k not set");
            sbwt::check_true(k <= KMER_MAX_LENGTH, "Maximum allowed k is " + std::to_string(KMER_MAX_LENGTH) + ". To increase the limit, recompile by first running cmake with the option `-DMAX_KMER_LENGTH=n`, where n is a number up to 255, and then running `make` again."); // 255 is max because of KMC
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

        if(from_index != ""){
            sbwt::check_true(colorfile == "", "Must not give both --from-index and --colorfile");
            sbwt::check_true(inputfile == "", "Must not give both --from-index and input sequences");
            sbwt::check_true(no_colors == false, "Must not give both --from-index and --load-dbg");
            sbwt::check_true(k == 0, "Must not give both --from-index and -k because k is defined in the index");
        }

        if(coloring_structure_type != "sdsl-fixed" && coloring_structure_type != "sdsl-hybrid" && coloring_structure_type != "roaring" && coloring_structure_type != "bitmagic"){
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
        ss << "Reverse complements = " << (reverse_complements ? "true" : "false") << "\n";
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

// A gzipped sequence reader with the possibility to reset the stream (required in color construction)
// This class exists because seekg(0) function in zstr::ifstream does not seem to work correctly.
class Gzip_Sequence_Reader_With_Reset{

    public:

        SeqIO::Reader<Buffered_ifstream<zstr::ifstream>>* reader;
        char* read_buf;
        int64_t read_buf_len;
        string filename;
        bool reverse_complements = false;

        Gzip_Sequence_Reader_With_Reset(string filename) : filename(filename){
            reader = new SeqIO::Reader<Buffered_ifstream<zstr::ifstream>>(filename);
            read_buf_len = 256;
            read_buf = (char*)malloc(read_buf_len);
        }

        int64_t get_next_read_to_buffer(){
            int64_t len = reader->get_next_read_to_buffer();

            // Copy the read from the internal buffer of the internal reader to our buffer
            if(len > read_buf_len){
                read_buf = (char*)realloc(read_buf, len);
                read_buf_len = len;
            }
            memcpy(read_buf, reader->read_buf, len);
            return len;
        }

        void enable_reverse_complements(){
            reverse_complements = true;
            reader->enable_reverse_complements();
        }

        void rewind_to_start(){
            delete(reader);
            reader = new SeqIO::Reader<Buffered_ifstream<zstr::ifstream>>(filename);
            if(reverse_complements) reader->enable_reverse_complements();
        }

        ~Gzip_Sequence_Reader_With_Reset(){
            delete reader;
            free(read_buf);
        }


};

// Returns filename of a new color file that has one color for each sequence
template<typename sequence_reader_t>
string generate_default_colorfile(sequence_reader_t& reader, bool reverse_complements){
    string colorfile = sbwt::get_temp_file_manager().create_filename();
    sbwt::Buffered_ofstream<> out(colorfile);
    stringstream ss;
    int64_t seq_id = 0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        ss.str(""); ss << seq_id << "\n";
        if(reverse_complements) ss << seq_id << "\n"; // Same color for the reverse complement pair
        out.write(ss.str().data(), ss.str().size());
        seq_id++;
    }
    return colorfile;
}


// Builds and serializes to disk
template<typename colorset_t> 
void build_coloring(plain_matrix_sbwt_t& dbg, const vector<int64_t>& color_assignment, const Build_Config& C){
    Coloring<colorset_t, typename colorset_t::view_t> coloring;
    if(C.input_format.gzipped){
        Coloring_Builder<colorset_t, typename colorset_t::view_t, Gzip_Sequence_Reader_With_Reset> cb; // Builder with gzipped input
        Gzip_Sequence_Reader_With_Reset reader(C.inputfile);
        if(C.reverse_complements) reader.enable_reverse_complements();
        cb.build_coloring(coloring, dbg, reader, color_assignment, C.memory_megas * (1 << 20), C.n_threads, C.colorset_sampling_distance);
    } else{
        Coloring_Builder<colorset_t, typename colorset_t::view_t, sbwt::SeqIO::Reader<Buffered_ifstream<std::ifstream>>> cb; // Builder without gzipped input
        sbwt::SeqIO::Reader<Buffered_ifstream<std::ifstream>> reader(C.inputfile);
        if(C.reverse_complements) reader.enable_reverse_complements();
        cb.build_coloring(coloring, dbg, reader, color_assignment, C.memory_megas * (1 << 20), C.n_threads, C.colorset_sampling_distance);        
    }
    sbwt::throwing_ofstream out(C.index_color_file, ios::binary);
    coloring.serialize(out.stream);
}

// Builds from existing index and serializes to disk
template<typename old_coloring_t, typename new_coloring_t> 
void build_from_index(plain_matrix_sbwt_t& dbg, const old_coloring_t& old_coloring, const Build_Config& C){

    // TODO: This makes a ton of unnecessary copies of things and has high peak RAM

    vector<typename new_coloring_t::colorset_type> new_colorsets;

    int64_t largest_color = 0;
    int64_t total_length = 0;

    for(int64_t i = 0; i < old_coloring.number_of_distinct_color_sets(); i++){
        vector<int64_t> set = old_coloring.get_color_set_as_vector_by_color_set_id(i);

        for(int64_t x : set) largest_color = max(largest_color, x);
        total_length += set.size();

        new_colorsets.push_back(set); // Re-encodes as new_coloring_t::color_set_t
    }

    typename new_coloring_t::colorset_storage_type new_storage(new_colorsets);
    new_coloring_t new_coloring(new_storage, 
                                old_coloring.get_node_id_to_colorset_id_structure(),
                                dbg, largest_color, total_length);

    throwing_ofstream colors_out(C.index_color_file);
    new_coloring.serialize(colors_out.stream);
    
    throwing_ofstream dbg_out(C.index_dbg_file);
    dbg.serialize(dbg_out.stream);

}

/*
    Coloring(const colorset_storage_type& sets,
             const Sparse_Uint_Array& node_id_to_color_set_id,
             const plain_matrix_sbwt_t& index,
             const int64_t largest_id,
             const int64_t total_color_set_length) 
             : sets(sets), node_id_to_color_set_id(node_id_to_color_set_id), index_ptr(&index), largest_color_id(largest_color_id), total_color_set_length(total_color_set_length){
    }

*/


// Creates a new file and returns the new filename.
// Transforms:
// 1
// 2
// 3
// Into:
// 1
// 1
// 2
// 2
// 3
// 3
string duplicate_each_color(const string& colorfile){
    Buffered_ifstream<> in(colorfile);

    string newfile = get_temp_file_manager().create_filename("",".txt");
    Buffered_ofstream<> out(newfile);

    char newline = '\n';
    string line;
    while(in.getline(line)){
        out.write(line.data(), line.size());
        out.write(&newline, 1);
        out.write(line.data(), line.size());
        out.write(&newline, 1);
    }

    return newfile;
}

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
        ("r,reverse-complements", "Also add reverse complements of the k-mers to the index.", cxxopts::value<bool>()->default_value("false"))
        ("temp-dir", "Directory for temporary files. This directory should have fast I/O operations and should have as much space as possible.", cxxopts::value<string>())
        ("m,mem-megas", "Number of megabytes allowed for external memory algorithms (must be at least 2048).", cxxopts::value<LL>()->default_value("2048"))
        ("t,n-threads", "Number of parallel exectuion threads. Default: 1", cxxopts::value<LL>()->default_value("1"))
        ("randomize-non-ACGT", "Replace non-ACGT letters with random nucleotides. If this option is not given, (k+1)-mers containing a non-ACGT character are deleted instead.", cxxopts::value<bool>()->default_value("false"))
        ("d,colorset-pointer-tradeoff", "This option controls a time-space tradeoff for storing and querying color sets. If given a value d, we store color set pointers only for every d nodes on every unitig. The higher the value of d, the smaller then index, but the slower the queries. The savings might be significant if the number of distinct color sets is small and the graph is large and has long unitigs.", cxxopts::value<LL>()->default_value("1"))
        ("no-colors", "Build only the de Bruijn graph without colors.", cxxopts::value<bool>()->default_value("false"))
        ("load-dbg", "If given, loads a precomputed de Bruijn graph from the index prefix. If this is given, the value of parameter -k is ignored because the order k is defined by the precomputed de Bruijn graph.", cxxopts::value<bool>()->default_value("false"))
        ("s,coloring-structure-type", "Type of coloring structure to build (\"sdsl-fixed\", \"sdsl-hybrid\", \"roaring\" or \"bitmagic\" ).", cxxopts::value<string>()->default_value("sdsl-hybrid"))
        ("from-index", "Take as input a pre-built Themisto index. Builds a new index in the format specified by --coloring-structure-type. This is currenlty implemented by decompressing the distinct color sets in memory before re-encoding them, so this might take a lot of RAM.", cxxopts::value<string>()->default_value(""))
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
    C.reverse_complements = opts["reverse-complements"].as<bool>();
    C.from_index = opts["from-index"].as<string>();

    if(C.verbose && C.silent) throw runtime_error("Can not give both --verbose and --silent");
    if(C.verbose) set_log_level(sbwt::LogLevel::MINOR);
    if(C.silent) set_log_level(sbwt::LogLevel::OFF);

    create_directory_if_does_not_exist(C.temp_dir);

    C.check_valid();
    sbwt::get_temp_file_manager().set_dir(C.temp_dir);

    write_log("Build configuration:\n" + C.to_string(), sbwt::LogLevel::MAJOR);
    write_log("Starting", sbwt::LogLevel::MAJOR);

    if(C.from_index != ""){
        // Transform existing index to new format

        sbwt::write_log("Loading de Bruijn Graph", sbwt::LogLevel::MAJOR);
        std::unique_ptr<sbwt::plain_matrix_sbwt_t> dbg_ptr = std::make_unique<sbwt::plain_matrix_sbwt_t>();
        dbg_ptr->load(C.index_dbg_file);

        sbwt::write_log("Loading coloring", sbwt::LogLevel::MAJOR);
        std::variant<Coloring<Color_Set, Color_Set_View>, Coloring<Roaring_Color_Set, Roaring_Color_Set>, Coloring<Bit_Magic_Color_Set, Bit_Magic_Color_Set>> old_coloring;
        load_coloring(C.index_color_file, *dbg_ptr, old_coloring);

        if(std::holds_alternative<Coloring<Color_Set, Color_Set_View>>(old_coloring))
            write_log("sdsl coloring structure loaded", LogLevel::MAJOR);
        if(std::holds_alternative<Coloring<Roaring_Color_Set, Roaring_Color_Set>>(old_coloring))
            write_log("roaring coloring structure loaded", LogLevel::MAJOR);
        if(std::holds_alternative<Coloring<Bit_Magic_Color_Set, Bit_Magic_Color_Set>>(old_coloring))
            write_log("BitMagic coloring structure loaded", LogLevel::MAJOR);

        auto visitor = [&](auto& old){
            if(C.coloring_structure_type == "sdsl-hybrid"){
                build_from_index<decltype(old), Coloring<Color_Set, Color_Set_View>>(*dbg_ptr, old, C);
            } else if(C.coloring_structure_type == "roaring"){
                build_from_index<decltype(old), Coloring<Roaring_Color_Set, Roaring_Color_Set>>(*dbg_ptr, old, C);
            } else if(C.coloring_structure_type == "bitmagic"){
                build_from_index<decltype(old), Coloring<Bit_Magic_Color_Set, Bit_Magic_Color_Set>>(*dbg_ptr, old, C);
            } else{
                throw std::runtime_error("Unkown coloring structure type: " + C.coloring_structure_type);
            }
        };

        std::visit(visitor, old_coloring);

        return 0;

    }

    // Deal with non-ACGT characters
    if(C.del_non_ACGT){
        // KMC takes care of this
    } else {
        write_log("Replacing non-ACGT characters with random nucleotides", LogLevel::MAJOR);
        C.inputfile = fix_alphabet(C.inputfile); // Turns the file into fasta format also
    }

    if(!C.no_colors){
        if(C.colorfile == ""){
            // Automatic colors
            sbwt::write_log("Assigning colors", sbwt::LogLevel::MAJOR);
            if(C.input_format.gzipped){
                sbwt::SeqIO::Reader<Buffered_ifstream<zstr::ifstream>> reader(C.inputfile);
                C.colorfile = generate_default_colorfile(reader, C.reverse_complements);
            } else{
                sbwt::SeqIO::Reader<Buffered_ifstream<std::ifstream>> reader(C.inputfile);
                C.colorfile = generate_default_colorfile(reader, C.reverse_complements);
            }
        } else{
            // User-defined colors
            if(C.reverse_complements){
                sbwt::write_log("Duplicating colors for reverse complements", sbwt::LogLevel::MAJOR);
                C.colorfile = duplicate_each_color(C.colorfile);
            }
        }
    }

    std::unique_ptr<sbwt::plain_matrix_sbwt_t> dbg_ptr;

    if(C.load_dbg){
        sbwt::write_log("Loading de Bruijn Graph", sbwt::LogLevel::MAJOR);
        dbg_ptr = std::make_unique<sbwt::plain_matrix_sbwt_t>();
        dbg_ptr->load(C.index_dbg_file);
    } else{
        sbwt::write_log("Building de Bruijn Graph", sbwt::LogLevel::MAJOR);

        vector<string> KMC_input_files = {C.inputfile};
        if(C.reverse_complements){
            write_log("Creating reverse complemented copy of " + C.inputfile + " to " + sbwt::get_temp_file_manager().get_dir(), LogLevel::MAJOR);
            sbwt::SeqIO::FileFormat fileformat = sbwt::SeqIO::figure_out_file_format(C.inputfile);;    
            if(fileformat.gzipped){
                KMC_input_files.push_back(sbwt::SeqIO::create_reverse_complement_files<
                    sbwt::SeqIO::Reader<sbwt::Buffered_ifstream<sbwt::zstr::ifstream>>,
                    sbwt::SeqIO::Writer<sbwt::Buffered_ofstream<sbwt::zstr::ofstream>>>({C.inputfile})[0]);
            } else{
                KMC_input_files.push_back(sbwt::SeqIO::create_reverse_complement_files<
                    sbwt::SeqIO::Reader<sbwt::Buffered_ifstream<std::ifstream>>,
                    sbwt::SeqIO::Writer<sbwt::Buffered_ofstream<std::ofstream>>>({C.inputfile})[0]);
            }
        }
        
        sbwt::plain_matrix_sbwt_t::BuildConfig sbwt_config;
        sbwt_config.build_streaming_support = true;
        sbwt_config.input_files = KMC_input_files;
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
        if(C.coloring_structure_type == "sdsl-hybrid"){
            build_coloring<Color_Set>(*dbg_ptr, color_assignment, C);
        } else if(C.coloring_structure_type == "roaring"){
            build_coloring<Roaring_Color_Set>(*dbg_ptr, color_assignment, C);
        } else if(C.coloring_structure_type == "bitmagic"){
            build_coloring<Bit_Magic_Color_Set>(*dbg_ptr, color_assignment, C);
        }
    } else{
        std::filesystem::remove(C.index_color_file); // There is an empty file so let's remove it
    }

    sbwt::write_log("Finished", sbwt::LogLevel::MAJOR);

    return 0;
}