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

class AlignerThread : public DispatcherConsumerCallback{

private:

    AlignerThread(const AlignerThread&); // No copying
    AlignerThread& operator=(const AlignerThread&); // No copying

    public:

        const plain_matrix_sbwt_t* SBWT; // Not owned by this class
        const Coloring* coloring; // Not owned by this class
        ParallelBaseWriter* out;
        bool reverse_complements;
        LL output_buffer_size;
        string output_buffer; // For printing
        LL k;

        vector<char> rc_buffer;
        

        AlignerThread(const plain_matrix_sbwt_t* SBWT, const Coloring* coloring, ParallelBaseWriter* out, bool reverse_complements, LL output_buffer_size){
            this->SBWT = SBWT;
            this->coloring = coloring;
            this->out = out;
            this->reverse_complements = reverse_complements;
            this->output_buffer_size = output_buffer_size;
            this->k = SBWT->get_k();
            rc_buffer.resize(1 << 10); // 1 kb. Will be resized if needed
        }


        virtual void callback(const char* S, LL S_size, int64_t string_id){
            if(S_size < k){
                write_log("Warning: query is shorter than k", LogLevel::MINOR);
                output_buffer += std::to_string(string_id) + "\n";
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

                // Todo: this needs memory allocation which can be bad for parallel performance?
                // Todo: faster int-to-string conversion by doing it manually
                output_buffer += std::to_string(string_id);
                for(color_t x : intersection.get_colors_as_vector()){
                    output_buffer += " " + std::to_string(x);
                }
                output_buffer += "\n";
            }

            // Flush buffer if needed
            if(output_buffer.size() > output_buffer_size){
                out->write(output_buffer);
                output_buffer.clear();
            }
            
        }

        virtual void finish(){
            if(output_buffer.size() > 0){
                out->write(output_buffer);
                output_buffer.clear();
            }
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