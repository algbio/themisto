#pragma once

#include "globals.hh"

Temp_File_Manager& get_temp_file_manager(){
    static Temp_File_Manager temp_file_manager; // Singleton
    return temp_file_manager;
}

long long cur_time_millis(){
	return (std::chrono::duration_cast< milliseconds >(system_clock::now().time_since_epoch())).count();
}

long long int program_start_millis = cur_time_millis();

double seconds_since_program_start(){
	return (cur_time_millis() - program_start_millis) / 1000.0;
}

string getTimeString(){
    std::time_t result = std::time(NULL);
    string time = std::asctime(std::localtime(&result));
    return time.substr(0,time.size() - 1); // Trim the trailing newline
}

bool logging_enabled = true;

void enable_logging(){
    logging_enabled = true;
}

void disable_logging(){
    logging_enabled = false;
}

std::mutex write_log_mutex;

void write_log(string message){
    std::lock_guard<std::mutex> lock(write_log_mutex);
    if(logging_enabled){
        std::streamsize default_precision = std::cout.precision();

        std::cerr << 
        std::setprecision(4) << std::fixed <<
        seconds_since_program_start() <<
        std::setprecision(default_precision) << 
        " " << getTimeString() << " " << message << std::endl;
    }
}

map<string,vector<string> > parse_args(int argc, char** argv){
    // Options are argumenta that start with "--". All non-options
    // that come after an option are parameters for that option
    map<string,vector<string> > M; // Option -> list of parameters
    string current_option = "";
    for(LL i = 1; i < argc; i++){
        string S = argv[i];
        if(S.size() >= 2  && S.substr(0,2) == "--"){
            current_option = S;
            M[current_option].resize(0); // Add empty vector for this option.
        } else{
            if(current_option == ""){
                cerr << "Error parsing command line parameters" << endl;
                exit(1);
            }
            M[current_option].push_back(S);
        }
    }
    return M;
}

string figure_out_file_format(string filename){
    for(LL i = (LL)filename.size()-1; i >= 0; i--){
        if(filename[i] == '.'){
            string end = filename.substr(i);
            
            if(end == ".fasta") return "fasta";
            if(end == ".fna") return "fasta";
            if(end == ".ffn") return "fasta";
            if(end == ".faa") return "fasta";
            if(end == ".frn") return "fasta";
            if(end == ".fa") return "fasta";

            if(end == ".fastq") return "fastq";
            if(end == ".fq") return "fastq";

            if(end == ".gz") return "gzip";

            throw(runtime_error("Unknown file format: " + filename));
        }
    }
    throw(runtime_error("Unknown file format: " + filename));
    return "unknown";
}

static constexpr char R_conv_tbl[] = { 'A', 'G' };
static constexpr char Y_conv_tbl[] = { 'C', 'T' };
static constexpr char K_conv_tbl[] = { 'G', 'T' };
static constexpr char M_conv_tbl[] = { 'A', 'C' };
static constexpr char S_conv_tbl[] = { 'C', 'G' };
static constexpr char W_conv_tbl[] = { 'A', 'T' };
static constexpr char B_conv_tbl[] = { 'C', 'G', 'T' };
static constexpr char D_conv_tbl[] = { 'A', 'G', 'T' };
static constexpr char H_conv_tbl[] = { 'A', 'C', 'T' };
static constexpr char V_conv_tbl[] = { 'A', 'C', 'G' };
static constexpr char N_conv_tbl[] = { 'A', 'C', 'G', 'T' };

char fix_char(char c){
	c = toupper(c);
    int rd = std::rand();
    
	switch (c) {
	case 'A':
		return c;
	case 'C':
		return c;
	case 'G':
		return c;
	case 'T':
		return c;
	case 'U':
		return 'T';
	case 'R':
		return R_conv_tbl[rd % 2];
	case 'Y':
		return Y_conv_tbl[rd % 2];
	case 'K':
		return K_conv_tbl[rd % 2];
	case 'M':
		return M_conv_tbl[rd % 2];
	case 'S':
		return S_conv_tbl[rd % 2];
	case 'W':
		return W_conv_tbl[rd % 2];
	case 'B':
		return B_conv_tbl[rd % 3];
	case 'D':
		return D_conv_tbl[rd % 3];
	case 'H':
		return H_conv_tbl[rd % 3];
	case 'V':
		return V_conv_tbl[rd % 3];
	default:
		return N_conv_tbl[rd % 4];
	}
}

// Returns number of chars replaced
LL fix_alphabet_of_string(string& S){
    LL chars_replaced = 0;
    for(LL i = 0; i < S.size(); i++){
        char c = S[i];
        char c_new = fix_char(c);
        if(c_new != c){
            S[i] = c_new;
            chars_replaced++;
        }
    }
    return chars_replaced;
}

// Makes a copy of the file and replaces bad characters. Returns the new filename
// The new file is in fasta format
std::string fix_alphabet(const std::string& input_file, const std::size_t bufsiz, const int mode){
    write_log("Making all characters upper case and replacing non-{A,C,G,T} characters with random characters from {A,C,G,T}");
    
    const std::string output_file = get_temp_file_manager().create_filename("seqs-");
    
    std::FILE* ip = std::fopen(input_file.c_str(), "rb");
    std::FILE* op = std::fopen(output_file.c_str(), "wb");
    
    char* ibuf = new char[bufsiz];
    char* obuf = new char[bufsiz];
    std::size_t sz;
    std::uint64_t replaced = 0;
    
    std::size_t j = 0;

    // FASTA
    if (mode == FASTA_MODE) {
        bool in_header = false;
        
        while ((sz = std::fread(ibuf, 1, bufsiz, ip))) {
            for (std::size_t i = 0; i < sz; ++i) {
                
                if (ibuf[i] == '>' || in_header) {
                    in_header = true;

                    if (ibuf[i] == '>') {
                        obuf[j++] = '>';
                    }
                    else if (ibuf[i] == '\n') {
                        in_header = false;
                        obuf[j++] = '\n';
                    }
                }
                
                else {
                    if (ibuf[i] == '\n') {
                        obuf[j++] = '\n';
                    }
                    else {
                        obuf[j] = fix_char(ibuf[i]);
                        if (obuf[j] != ibuf[i]) {
                            ++replaced;
                        }
                        ++j;
                    }
                }
            }
            std::fwrite(obuf, 1, j, op);
            j = 0;
        }
    }
    // FASTQ
    else {
        enum FASTQ_line { seqname, seq, plus, qual };
        FASTQ_line fl = seqname;

        while ((sz = std::fread(ibuf, 1, bufsiz, ip))) {
            
            for (std::size_t i = 0; i < sz; ++i) {

                switch(fl) {
                    
                case FASTQ_line::seqname : {
                    if (ibuf[i] == '@') {
                        obuf[j++] = '>';
                    }
                    else if (ibuf[i] == '\n') {
                        fl = FASTQ_line::seq;
                        obuf[j++] = '\n';
                    }
                }
                    break;

                case FASTQ_line::seq : {
                    if (ibuf[i] == '\n') {
                        fl = FASTQ_line::plus;
                        obuf[j++] = '\n';
                    }
                    else {
                        obuf[j] = fix_char(ibuf[i]);
                        if (obuf[j] != ibuf[i]) {
                            ++replaced;
                        }
                        ++j;
                    }
                }
                    break;

                case FASTQ_line::plus : {
                    if (ibuf[i] == '\n') {
                        fl = FASTQ_line::qual;
                    }
                }
                    break;

                case FASTQ_line::qual : {
                    if (ibuf[i] == '\n') {
                        fl = FASTQ_line::seqname;
                    }
                }
                    break;
                    
                default :
                    break;
                }
            }
            std::fwrite(obuf, 1, j, op);
            j = 0;
        }
    }

    delete[] ibuf;
    delete[] obuf;
    
    std::fclose(ip);
    std::fclose(op);
    
    write_log("Replaced " + to_string(replaced) + " characters");
    
    return output_file;
}

// We need this function because the standard library stoll function accepts all kinds of crap,
// such as "123aasfjhk" and "1 2 3 4" as a number. This function check that the string is a valid
// number and returns that number, or throws an error otherwise.
LL string_to_integer_safe(const string& S){

    // Figure out leading and trailing whitespace
    LL pos_of_first_digit = 1e18;
    LL pos_of_last_digit = -1;
    for(LL i = 0; i < (LL)S.size(); i++){
        if(!std::isdigit(S[i]) && !std::isspace(S[i]))
            throw std::runtime_error("Error parsing color file: could not parse integer: " + S);
        if(std::isdigit(S[i])){
            pos_of_first_digit = min(pos_of_first_digit, i);
            pos_of_last_digit = max(pos_of_last_digit, i);
        }
    }
    if(pos_of_last_digit == -1) throw std::runtime_error("Error parsing color file: could not parse integer: " + S); // No digits found

    // Check that there are no internal spaces
    for(LL i = pos_of_first_digit; i <= pos_of_last_digit; i++)
        if(!std::isdigit(S[i])) throw std::runtime_error("Error parsing color file: could not parse integer: " + S); // Internat whitespace

    // Checks ok, convert to integer
    return stoll(S);
}

vector<LL> read_colorfile(string filename){
    vector<LL> seq_to_color;
    throwing_ifstream colors_in(filename);
    string line;
    while(colors_in.getline(line)){
        seq_to_color.push_back(string_to_integer_safe(line));
    }
    return seq_to_color;
}

// Returns new inputfile and new colorfile
pair<string,string> split_all_seqs_at_non_ACGT(string inputfile, string inputfile_format, string colorfile){
    string new_colorfile = get_temp_file_manager().create_filename();
    string new_seqfile = get_temp_file_manager().create_filename();

    throwing_ofstream colors_out(new_colorfile);

    Sequence_Reader fr(inputfile, inputfile_format == "fasta" ? FASTA_MODE : FASTQ_MODE);
    throwing_ofstream sequences_out(new_seqfile);
    vector<LL> colors = read_colorfile(colorfile);
    LL seq_id = 0;
    while(!fr.done()){

        // Read a sequence and its color
        string seq = fr.get_next_query_stream().get_all();

        // Chop the sequence into pieces that have only ACGT characters
        seq += '$'; // Trick to avoid having a special case for the last sequence
        string new_seq;
        for(char c : seq){
            assert((c >= 'A' && c <= 'Z') || c == '$');
            if(c == 'A' || c == 'C' || c == 'G' || c == 'T')
                new_seq += c;
            else{
                if(new_seq.size() > 0){
                    sequences_out << ">\n" << new_seq << "\n";
                    colors_out << colors[seq_id] << "\n";
                    new_seq = "";
                }
            }
        }
        seq_id++;
    }
    return {new_seqfile, new_colorfile};
}



void sigint_handler(int sig) {
    cerr << "caught signal: " << sig << endl;
    cerr << "Cleaning up temporary files" << endl;
    get_temp_file_manager().delete_all_files();
    exit(1);
}

void sigabrt_handler(int sig) {
    cerr << "caught signal: " << sig << endl;
    cerr << "Cleaning up temporary files" << endl;
    get_temp_file_manager().delete_all_files();
    cerr << "Aborting" << endl;
    exit(1);
}

auto sigint_register_return_value = signal(SIGINT, sigint_handler); // Set the SIGINT handler
auto sigabrt_register_return_value = signal(SIGABRT, sigabrt_handler); // Set the SIGABRT handler

vector<string> get_first_and_last_kmers(string fastafile, LL k){
    // todo: this is pretty expensive because this has to read the whole reference data
    Sequence_Reader fr(fastafile, FASTA_MODE);
    vector<string> result;
    while(!fr.done()){
        string ref = fr.get_next_query_stream().get_all();
        if(ref.size() >= k){
            result.push_back(ref.substr(0,k));
            result.push_back(ref.substr(ref.size()-k,k));
        }
    }
    return result;
}

string get_rc(string S){
    std::reverse(S.begin(), S.end());
    for(char& c : S){
        if(c == 'A') c = 'T';
        else if(c == 'C') c = 'G';
        else if(c == 'G') c = 'C';
        else if(c == 'T') c = 'A';
    }
    return S;
}

// true if S is colexicographically-smaller than T
bool colex_compare(const string& S, const string& T){
    LL i = 0;
    while(true){
        if(i == S.size() || i == T.size()){
            // One of the strings is a suffix of the other. Return the shorter.
            if(S.size() < T.size()) return true;
            else return false;
        }
        if(S[S.size()-1-i] < T[T.size()-1-i]) return true;
        if(S[S.size()-1-i] > T[T.size()-1-i]) return false;
        i++;
    }
}

bool colex_compare_cstrings(const char* x, const char* y){
    LL nx = strlen(x);
    LL ny = strlen(y);
    for(LL i = 0; i < min(nx,ny); i++){
        if(x[nx-1-i] < y[ny-1-i]) return true;
        if(x[nx-1-i] > y[ny-1-i]) return false;
    }

    // All no mismatches -> the shorter string is smaller
    return nx < ny;
};

bool lex_compare(const string& S, const string& T){
    return S < T;
};

bool lex_compare_cstrings(const char* x, const char* y){
    return strcmp(x,y) < 0;
};

// Split by whitespace
vector<string> split(string text){
    std::istringstream iss(text);
    std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                 std::istream_iterator<std::string>());
    return results;
}

// Split by delimiter
vector<string> split(string text, char delimiter){
    assert(text.size() != 0); // If called with empty string we probably have a bug
    vector<LL> I; // Delimiter indices
    I.push_back(-1);
    for(LL i = 0; i < text.size(); i++){
        if(text[i] == delimiter){
            I.push_back(i);
        }
    }
    I.push_back(text.size());
    vector<string> tokens;
    for(LL i = 0; i < I.size()-1; i++){
        LL len = I[i+1] - I[i] + 1 - 2;
        tokens.push_back(text.substr(I[i]+1, len));
    }
    
    return tokens;
}

vector<string> split(const char* text, char delimiter){
    return split(string(text), delimiter);
}

void create_directory_if_does_not_exist(string path){
    std::filesystem::create_directory(path);
}

// https://stackoverflow.com/questions/18100097/portable-way-to-check-if-directory-exists-windows-linux-c
void check_dir_exists(string path){
    struct stat info;    
    if( stat( path.c_str(), &info ) != 0 ){
        cerr << "Error: can not access directory " << path << endl;
        exit(1);
    }
    else if( info.st_mode & S_IFDIR ){
        // All good
    }    
    else{
        cerr << "Error: is not a directory: " << path << endl;
        exit(1);
    }
}

void check_readable(string filename){
    throwing_ifstream F(filename); // Throws on failure
}

// Also clears the file
void check_writable(string filename){
    throwing_ofstream F(filename, std::ofstream::out | std::ofstream::app); // Throws on failure
}

vector<string> get_all_lines(string infile){
    vector<string> lines;
    string line;
    throwing_ifstream in(infile);
    while(in.getline(line)){
        lines.push_back(line);
    }
    return lines;
}

vector<char> read_binary_file(string infile){
    throwing_ifstream file(infile, std::ios::binary | std::ios::ate);
    std::streamsize size = file.stream.tellg();
    file.stream.seekg(0, std::ios::beg);

    std::vector<char> buffer(size);
    if (file.read(buffer.data(), size)){
        return buffer;
    } else{
        cerr << "Error reading file: " << infile << endl;
        assert(false);
    }
}

bool files_are_equal(const std::string& p1, const std::string& p2) {
  //https://stackoverflow.com/questions/6163611/compare-two-files/6163627
    throwing_ifstream f1(p1, std::ifstream::binary|std::ifstream::ate);
    throwing_ifstream f2(p2, std::ifstream::binary|std::ifstream::ate);

    if (f1.stream.tellg() != f2.stream.tellg()) {
      return false; //size mismatch
    }

    //seek back to beginning and use std::equal to compare contents
    f1.stream.seekg(0, std::ifstream::beg);
    f2.stream.seekg(0, std::ifstream::beg);
    return std::equal(std::istreambuf_iterator<char>(f1.stream.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.stream.rdbuf()));
}

void check_true(bool condition, string error_message){
    if(!condition){
        throw std::runtime_error(error_message);
    }
}

// Returns filename of a new color file that has one color for each sequence
// Input format is either "fasta" or "fastq"
string generate_default_colorfile(string inputfile, string file_format){
    string colorfile = get_temp_file_manager().create_filename();
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