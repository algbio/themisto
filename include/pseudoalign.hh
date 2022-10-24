#pragma once

#include <string>
#include "sbwt/SBWT.hh"
#include "new_coloring.hh"
#include "SeqIO.hh"
#include "variants.hh"

using namespace std;
using namespace sbwt;

typedef int64_t color_t;

// ASSUMES x >= -1
// Buffer must have enough space to fit the ascii representation of integer x
// Returns the length of the string written to the buffer.
int64_t fast_int_to_string(int64_t x, char* buffer);

// returns a vector where element i is the ref ids aligned with query i
vector<set<LL> > parse_pseudoalignment_output_format_from_disk(string filename);


template<class coloring_t>
class AlignerThread : public DispatcherConsumerCallback{

private:

    AlignerThread(const AlignerThread&); // No copying
    AlignerThread& operator=(const AlignerThread&); // No copying

    public:

        const plain_matrix_sbwt_t* SBWT; // Not owned by this class
        const coloring_t* coloring; // Not owned by this class
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

        // Buffer for storing color set ids
        vector<int64_t> color_set_id_buffer;
        vector<int64_t> rc_color_set_id_buffer;

        AlignerThread(const plain_matrix_sbwt_t* SBWT, const coloring_t* coloring, ParallelBaseWriter* out, bool reverse_complements, LL output_buffer_capacity){
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

        // -1 if node is not found at all.
        void push_color_set_ids_to_buffer(const vector<int64_t>& colex_ranks, vector<int64_t>& buffer){

            // First pass: get all k-mer kmers and the last k-mer
            for(int64_t v : colex_ranks){
                if(v == -1) buffer.push_back(-1); // k-mer not found
                else{
                    if(coloring->is_core_kmer(v) || v == colex_ranks.back()){
                        int64_t id = coloring->get_color_set_id(v);
                        buffer.push_back(id);
                    } else{
                        buffer.push_back(-2); // To be filled in the second pass
                    }
                }
            }

            // Second pass: fill in the rest
            for(int64_t i = (int64_t)colex_ranks.size()-2; i >= 0; i--){ // -2: skip the last one
                if(buffer[i] == -2){
                    if(buffer[i+1] == -1){
                        // Can't copy from the next k-mer because it's not found
                        buffer[i] = coloring->get_color_set_id(colex_ranks[i]);
                    }
                    else {
                        // Can copy from the next k-mer
                        buffer[i] = buffer[i+1];
                    }
                }
            }
        }

        // Returns the color set
        vector<int64_t> do_intersections_on_color_id_buffers_with_reverse_complements(){
            LL n_kmers = color_set_id_buffer.size();

            bool first_nonempty_union_found = false;
            typename coloring_t::colorset_type result;
            for(LL i = 0; i < n_kmers; i++){
                if(i > 0
                && (color_set_id_buffer[i] == color_set_id_buffer[i-1])
                && (rc_color_set_id_buffer[n_kmers-1-i] == rc_color_set_id_buffer[n_kmers-1-i+1])){
                    continue; // This pair of color set ids was already intersected in the previous iteration
                }

                int64_t fw_id = color_set_id_buffer[i];
                int64_t rc_id = rc_color_set_id_buffer[n_kmers-1-i];

                typename coloring_t::colorset_type cs; // For union. TODO: Reuse a union buffer
                if(fw_id == -1 && rc_id == -1) continue; // Neither direction is found
                else if(fw_id == -1 && rc_id >= 0) cs = coloring->get_color_set_by_color_set_id(rc_id);
                else if(fw_id >= 0 && rc_id == -1) cs = coloring->get_color_set_by_color_set_id(fw_id);
                else if(fw_id >= 0 && rc_id >= 0){
                    // Take union of forward and reverse complement
                    cs = coloring->get_color_set_by_color_set_id(fw_id);
                    cs.do_union(coloring->get_color_set_by_color_set_id(rc_id));
                }

                if(!cs.empty()){
                    if(!first_nonempty_union_found){
                        result = cs; // This is the first nonempty union
                        first_nonempty_union_found = true;
                    }
                    else
                        result.intersection(cs); // Intersection
                }
            }
            return result.get_colors_as_vector();
        }

        // Returns the color se
        vector<int64_t> do_intersections_on_color_id_buffers_without_reverse_complements(){
            LL n_kmers = color_set_id_buffer.size();

            bool first_nonempty_color_set_found = false;
            typename coloring_t::colorset_type result;
            for(LL i = 0; i < n_kmers; i++){
                if(i > 0  && (color_set_id_buffer[i] == color_set_id_buffer[i-1])){
                    continue; // This color set was already intersected in the previous iteration
                }
                if(color_set_id_buffer[i] == -1) continue; // k-mer not found

                const typename coloring_t::colorset_type& cs = coloring->get_color_set_by_color_set_id(color_set_id_buffer[i]);
                if(cs.size() > 0){
                    if(!first_nonempty_color_set_found){
                        result = cs; // This is the first nonempty color set
                        first_nonempty_color_set_found = true;
                    }
                    else
                        result.intersection(cs); // Intersection
                }
            }
            return result.get_colors_as_vector();
        }

        virtual void callback(const char* S, LL S_size, int64_t string_id){
            char string_to_int_buffer[32]; // Enough space for a 64-bit integer in ascii
            char newline = '\n';
            char space = ' ';

            color_set_id_buffer.resize(0);
            rc_color_set_id_buffer.resize(0);
            // Clearing the buffers like this might look bad for performance at first glance because
            // we will then need to allocate new space for new elements that will be pushed to the buffers.
            // But in fact it's ok because resize is not supposed to affect the internal capacity of the vector.
            // cppreference.com says:
            //     "Vector capacity is never reduced when resizing to smaller size because that would
            //      invalidate all iterators, rather than only the ones that would be invalidated by the
            //      equivalent sequence of pop_back() calls."

            if(S_size < k){
                write_log("Warning: query is shorter than k", LogLevel::MINOR);
                int64_t len = fast_int_to_string(string_id, string_to_int_buffer);
                add_to_output(string_to_int_buffer, len);
                add_to_output(&newline, 1);
            }
            else{
                vector<int64_t> colex_ranks = SBWT->streaming_search(S, S_size); // TODO: version that pushes to existing buffer?
                push_color_set_ids_to_buffer(colex_ranks, color_set_id_buffer);
                vector<int64_t> rc_colex_ranks;
                if(reverse_complements){
                    while(S_size > rc_buffer.size()){ // Make sure buffer is long enough. TODO: cleaner code with resize()?
                        rc_buffer.resize(rc_buffer.size()*2);
                    }
                    memcpy(rc_buffer.data(), S, S_size);
                    reverse_complement_c_string(rc_buffer.data(), S_size); // There is no null at the end but that is ok
                    rc_colex_ranks = SBWT->streaming_search(rc_buffer.data(), S_size);
                    push_color_set_ids_to_buffer(rc_colex_ranks, rc_color_set_id_buffer);
                }

                vector<int64_t> intersection;
                if(reverse_complements) intersection = do_intersections_on_color_id_buffers_with_reverse_complements();
                else intersection = do_intersections_on_color_id_buffers_without_reverse_complements();

                int64_t len = fast_int_to_string(string_id, string_to_int_buffer);
                add_to_output(string_to_int_buffer, len);
                for(color_t x : intersection){
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

template<typename coloring_t, typename sequence_reader_t>
void pseudoalign(const plain_matrix_sbwt_t& SBWT, const coloring_t& coloring, int64_t n_threads, sequence_reader_t& reader, std::string outfile, bool reverse_complements, int64_t buffer_size, bool gzipped, bool sorted_output){

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
        AlignerThread<coloring_t>* T = new AlignerThread<coloring_t>(&SBWT, &coloring, out, reverse_complements, buffer_size);
        threads.push_back(T);
    }

    run_dispatcher(threads, reader, buffer_size);

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
