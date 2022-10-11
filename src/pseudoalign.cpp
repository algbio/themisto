#include <string>
#include <sstream>
#include "pseudoalign.hh"
#include "WorkDispatcher.hh"

using namespace std;
using namespace sbwt;

LL intersect_buffers(vector<color_t>& buf1, LL buf1_len, vector<color_t>& buf2, LL buf2_len){

    LL i = 0, j = 0, k = 0;
    while(i < buf1_len && j < buf2_len){
        if(buf1[i] < buf2[j]) i++;
        else if(buf1[i] > buf2[j]) j++;
        else{
            buf1[k] = buf1[i];
            i++; j++; k++;
        }
    }
    return k;

}

LL union_buffers(vector<color_t>& buf1, LL buf1_len, vector<color_t>& buf2, LL buf2_len, vector<color_t>& result_buf){

    auto end = std::set_union(
                    buf1.begin(), buf1.begin() + buf1_len, 
                    buf2.begin(), buf2.begin() + buf2_len, 
                    result_buf.begin()
    );
    return end - result_buf.begin();
}

template <typename T>
string vec_to_string(const vector<T>& v){
    stringstream ss;
    for(T x : v) ss << x << " ";
    return ss.str();
}

vector<set<LL> > parse_pseudoalignment_output_format_from_disk(string filename){
    vector<pair<LL, set<LL>>> results; // Pairs (query id, color set)
    check_readable(filename);
    throwing_ifstream input(filename);
    string line;
    
    LL line_number = 0;
    while(input.getline(line)){
        vector<LL> tokens = parse_tokens<LL>(line);
        assert(tokens.size() >= 1);
        LL query_id = tokens[0];
        set<LL> alignments;
        for(LL i = 1; i < tokens.size(); i++){
            alignments.insert(tokens[i]);
        }
        results.push_back({query_id, alignments});
        line_number++;
    }

    sort(results.begin(), results.end());
    vector<set<LL> > just_color_sets;
    for(auto X : results) just_color_sets.push_back(X.second);
    return just_color_sets;
}

// ASSUMES x >= -1
// Buffer must have enough space to fit the ascii representation of integer x
// Returns the length of the string written to the buffer.
int64_t fast_int_to_string(int64_t x, char* buffer){
    // Fast manual integer-to-string conversion
    
    LL i = 0;
    // Write the digits in reverse order (reversed back at the end)
    if(x == -1){
        buffer[0] = '1';
        buffer[1] = '-';
        i = 2;
    } else if(x == 0){
        buffer[0] = '0';
        i = 1;
    } else{
        while(x > 0){
            buffer[i++] = '0' + (x % 10);
            x /= 10;
        }
    }
    std::reverse(buffer, buffer + i);
    buffer[i] = '\0';
    return i;
}

class AlignerThread : public DispatcherConsumerCallback{

private:

    AlignerThread(const AlignerThread&); // No copying
    AlignerThread& operator=(const AlignerThread&); // No copying

    public:

        const plain_matrix_sbwt_t* SBWT; // Not owned by this class
        const Coloring* coloring; // Not owned by this class
        ParallelBaseWriter* out;
        bool reverse_complements;
        LL output_buffer_capacity;
        LL output_buffer_size;
        LL k;

        // Buffer for reverse-complementing strings
        vector<char> rc_buffer;

        // Buffer for printing. We want to have a local buffer for each thread to avoid having to call the 
        // parallel writer so often to avoid locking the writer from other threads.
        vector<char> output_buffer;

        AlignerThread(const plain_matrix_sbwt_t* SBWT, const Coloring* coloring, ParallelBaseWriter* out, bool reverse_complements, LL output_buffer_capacity){
            this->SBWT = SBWT;
            this->coloring = coloring;
            this->out = out;
            this->reverse_complements = reverse_complements;
            this->output_buffer_capacity = output_buffer_capacity;
            this->k = SBWT->get_k();
            rc_buffer.resize(1 << 10); // 1 kb. Will be resized if needed
            output_buffer.resize(output_buffer_capacity);
            output_buffer_size = 0;
        }

        void add_to_output(char* data, int64_t data_length){
            if(data_length > output_buffer_capacity) throw std::runtime_error("Bug: output buffer too small");
            
            if(output_buffer_size + data_length > output_buffer_capacity){
                // Flush the buffer
                out->write(output_buffer.data(), output_buffer_size);
                output_buffer_size = 0;
            }
            memcpy(output_buffer.data() + output_buffer_size, data, data_length); // Move data to buffer
            output_buffer_size += data_length;
        }


        virtual void callback(const char* S, LL S_size, int64_t string_id){
            char string_to_int_buffer[32]; // Enough space for a 64-bit integer in ascii
            char newline = '\n';
            char space = ' ';

            if(S_size < k){
                write_log("Warning: query is shorter than k", LogLevel::MINOR);
                int64_t len = fast_int_to_string(string_id, string_to_int_buffer);
                add_to_output(string_to_int_buffer, len);
                add_to_output(&newline, 1);
            }
            else{
                vector<int64_t> colex_ranks = SBWT->streaming_search(S, S_size);
                vector<int64_t> rc_colex_ranks;
                if(reverse_complements){

                    while(S_size > rc_buffer.size()){ // Make sure buffer is long enough
                        rc_buffer.resize(rc_buffer.size()*2);
                    }
                    memcpy(rc_buffer.data(), S, S_size);
                    reverse_complement_c_string(rc_buffer.data(), S_size); // There is no null at the end but that is ok
                    rc_colex_ranks = SBWT->streaming_search(rc_buffer.data(), S_size);
                }
                LL n_kmers = colex_ranks.size();
                Color_Set intersection;
                bool first_nonempty_found = false;
                for(LL i = 0; i < n_kmers; i++){
                    if(colex_ranks[i] >= 0 || (reverse_complements && rc_colex_ranks[n_kmers-1-i] >= 0)){ // k-mer found
                        Color_Set cs; // Empty
                        if(colex_ranks[i] >= 0) cs = coloring->get_color_set(colex_ranks[i]);
                        if(reverse_complements && rc_colex_ranks[n_kmers-1-i] >= 0){
                            Color_Set cs_rc = coloring->get_color_set(rc_colex_ranks[n_kmers-1-i]);
                            cs = cs.do_union(cs_rc);
                        }
                        if(cs.size() > 0){
                            if(!first_nonempty_found){
                                intersection = cs;
                                first_nonempty_found = true;
                            } else{
                                intersection = intersection.intersection(cs);
                            }
                        }
                    }
                }

                int64_t len = fast_int_to_string(string_id, string_to_int_buffer);
                add_to_output(string_to_int_buffer, len);
                for(color_t x : intersection.get_colors_as_vector()){
                    len = fast_int_to_string(x, string_to_int_buffer);
                    add_to_output(&space, 1);
                    add_to_output(string_to_int_buffer, len);
                }
                add_to_output(&newline, 1);
            }            
        }

        virtual void finish(){
            out->write(output_buffer.data(), output_buffer_size);
            output_buffer_size = 0;
        }
};

template<typename instream_t, typename outstream_t>
void sort_parallel_output_file(instream_t& instream, outstream_t& outstream){
    set<pair<LL, string> > Q; // Priority queue with pairs (priority, content)
    string line;
    vector<string> tokens;
    LL current_query_id = 0;
    
    while(getline(instream,line)){
        stringstream ss(line);
        LL priority; ss >> priority;
        Q.insert({priority, line + "\n"});
        while(Q.begin()->first == current_query_id){
            outstream << Q.begin()->second;
            Q.erase(Q.begin());
            current_query_id++;
        }
    }
    assert(Q.empty());
}

void pseudoalign(const plain_matrix_sbwt_t& SBWT, const Coloring& coloring, int64_t n_threads, std::string inputfile, std::string outfile, bool reverse_complements, int64_t buffer_size, bool gzipped, bool sorted_output){

    ParallelBaseWriter* out = nullptr;
    if (!outfile.empty()) {
        if (gzipped)
            out = new ParallelGzipWriter(outfile);
        else
            out = new ParallelOutputWriter(outfile);
    } else {
        if (gzipped)
            out = new ParallelGzipWriter(cout);
        else
            out = new ParallelOutputWriter(cout);
    }

    vector<DispatcherConsumerCallback*> threads;
    for (LL i = 0; i < n_threads; i++) {
        AlignerThread* T = new AlignerThread(&SBWT, &coloring, out, reverse_complements, buffer_size);
        threads.push_back(T);
    }

    sbwt::SeqIO::Reader<> sr(inputfile);
    run_dispatcher(threads, sr, buffer_size);

    // Clean up
    for (DispatcherConsumerCallback* t : threads) delete t;
    out->flush();
    delete out;

    if (sorted_output) {
        write_log("Sorting output file", LogLevel::MAJOR);
        string tempfile =
            get_temp_file_manager().create_filename("results_temp");
        if (gzipped) {
            zstr::ifstream instream(outfile);
            zstr::ofstream outstream(tempfile);
            sort_parallel_output_file(instream, outstream);
        } else {
            throwing_ifstream instream(outfile);
            throwing_ofstream outstream(tempfile);
            sort_parallel_output_file(instream.stream, outstream.stream);
        }
        std::filesystem::rename(tempfile, outfile);
    }
}