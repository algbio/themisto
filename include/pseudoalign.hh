#pragma once

#include <string>
#include "sbwt/SBWT.hh"
#include "coloring/Coloring.hh"
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
vector<vector<int64_t> >  parse_pseudoalignment_output_format_from_disk(string filename);

template<class coloring_t>
class Pseudoaligner_Base{

    // This class contains code that is common to both the intersection and threshold pseudualigners

public:

        const plain_matrix_sbwt_t* SBWT; // Not owned by this class
        const coloring_t* coloring; // Not owned by this class
        ParallelBaseWriter* out;
        bool reverse_complements;
        int64_t k;

        // Buffer for reverse-complementing strings
        vector<char> rc_buffer;

        // Buffer for printing. We want to have a local buffer for each thread to avoid having to call the
        // parallel writer so often to avoid locking the writer from other threads.
        vector<char> output_buffer;
        int64_t output_buffer_flush_threshold;

        // Buffer for storing color set ids
        vector<int64_t> color_set_id_buffer;
        vector<int64_t> rc_color_set_id_buffer;

        Pseudoaligner_Base(const plain_matrix_sbwt_t* SBWT, const coloring_t* coloring, ParallelBaseWriter* out, bool reverse_complements, int64_t output_buffer_capacity){
            this->SBWT = SBWT;
            this->coloring = coloring;
            this->out = out;
            this->reverse_complements = reverse_complements;
            this->output_buffer_flush_threshold = output_buffer_capacity;
            this->k = SBWT->get_k();
            rc_buffer.resize(1 << 10); // 1 kb. Will be resized if needed
            output_buffer.reserve(output_buffer_capacity);
        }

        // Flushes at the next newline when output_buffer_flush_threshold is exceeded
        void add_to_output(char* data, int64_t data_length){
            for(int64_t i = 0; i < data_length; i++){
                output_buffer.push_back(data[i]);
                if(output_buffer.size() > output_buffer_flush_threshold && data[i] == '\n'){
                    out->write(output_buffer.data(), output_buffer.size());
                    output_buffer.clear(); // Let's hope this keeps the reserved capacity of the vector intact
                }
            }
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

        vector<int64_t> get_rc_colex_ranks(const char* S, int64_t S_size){
            while(S_size > rc_buffer.size()){
                rc_buffer.resize(rc_buffer.size()*2);
            }
            memcpy(rc_buffer.data(), S, S_size);
            reverse_complement_c_string(rc_buffer.data(), S_size); // There is no null at the end but that is ok
            return SBWT->streaming_search(rc_buffer.data(), S_size);
        }

};

template<class coloring_t>
class ThresholdPseudoaligner : public DispatcherConsumerCallback, Pseudoaligner_Base<coloring_t>{

private:

typedef Pseudoaligner_Base<coloring_t> Base;
double count_threshold; // Fraction of k-mers that need to be found to report pseudoalignment to a color
bool ignore_unknown_kmers = false; // Ignore k-mers that do not exist in the de Bruijn graph or have no colors

// State used during callback
vector<int64_t> counts; // counts[i] = number of occurrences of color i
vector<int64_t> nonzero_count_indices; // Indices in this->counts that have a non-zero value

public:

    ThresholdPseudoaligner(const plain_matrix_sbwt_t* SBWT, const coloring_t* coloring, ParallelBaseWriter* out, bool reverse_complements, int64_t output_buffer_capacity, double count_threshold, bool ignore_unknown) : Pseudoaligner_Base<coloring_t>(SBWT, coloring, out, reverse_complements, output_buffer_capacity), count_threshold(count_threshold), ignore_unknown_kmers(ignore_unknown){
        counts.resize(coloring->largest_color() + 1); // Initializes counts to zeroes
    }

    virtual void callback(const char* S, int64_t S_size, int64_t string_id, std::array<uint8_t, 8> metadata){
        (void) metadata; // Ignored

        char string_to_int_buffer[32]; // Enough space for a 64-bit integer in ascii
        char newline = '\n';
        char space = ' ';

        Base::color_set_id_buffer.resize(0);
        Base::rc_color_set_id_buffer.resize(0);

        if(S_size < Base::k){
            write_log("Warning: query is shorter than k", LogLevel::MINOR);
            int64_t len = fast_int_to_string(string_id, string_to_int_buffer);
            Base::add_to_output(string_to_int_buffer, len);
            Base::add_to_output(&newline, 1);
        } else{
            vector<int64_t> colex_ranks = Base::SBWT->streaming_search(S, S_size); // TODO: version that pushes to existing buffer?
            Base::push_color_set_ids_to_buffer(colex_ranks, Base::color_set_id_buffer);
            if(Base::reverse_complements){
                Base::push_color_set_ids_to_buffer(Base::get_rc_colex_ranks(S, S_size), Base::rc_color_set_id_buffer);
            }

            typename coloring_t::colorset_type fw_set;
            int64_t n_kmers = S_size - Base::k  + 1;
            int64_t n_kmers_with_at_least_1_color = 0;
            int64_t run_length = 0; // Number of consecutive identical color sets 
            for(int64_t kmer_idx = 0; kmer_idx < n_kmers; kmer_idx++){
                run_length++;

                bool last = (kmer_idx == n_kmers - 1);
                bool fw_different = !last && Base::color_set_id_buffer[kmer_idx] != Base::color_set_id_buffer[kmer_idx+1];
                bool rc_different = !last && Base::reverse_complements && Base::rc_color_set_id_buffer[n_kmers-1-kmer_idx] != Base::rc_color_set_id_buffer[n_kmers-1-kmer_idx-1];
                bool end_of_run = (last || fw_different || rc_different);

                if(end_of_run){

                    // Retrieve forward color set
                    if(Base::color_set_id_buffer[kmer_idx] == -1) fw_set = (typename coloring_t::colorset_type){}; // Empty
                    else fw_set = Base::coloring->get_color_set_by_color_set_id(Base::color_set_id_buffer[kmer_idx]);
                    
                    // Retrieve reverse complement color set
                    if(Base::reverse_complements && Base::rc_color_set_id_buffer[n_kmers - 1  - kmer_idx] != -1){
                        typename coloring_t::colorset_type::view_t rc_set_view = Base::coloring->get_color_set_by_color_set_id(Base::rc_color_set_id_buffer[n_kmers - 1 - kmer_idx]);
                        fw_set.do_union(rc_set_view);
                    }

                    bool has_at_least_one_color =false;
                    // Add the run length to the counts
                    for(int64_t color : fw_set.get_colors_as_vector()){
                        has_at_least_one_color = true;
                        if(counts[color] == 0) nonzero_count_indices.push_back(color);
                        counts[color] += run_length;
                    }

                    n_kmers_with_at_least_1_color += has_at_least_one_color * run_length;

                    run_length = 0; // Reset the run
                }
            }

            // Print the colors of all counters that are above threshold
            int64_t len = fast_int_to_string(string_id, string_to_int_buffer);
            Base::add_to_output(string_to_int_buffer, len); // String id
            for(int64_t color : nonzero_count_indices){
                int64_t count = counts[color];
                int64_t effective_kmers = ignore_unknown_kmers ? n_kmers_with_at_least_1_color : n_kmers;
                if(count >= effective_kmers * count_threshold){
                    // Report color
                    len = fast_int_to_string(color, string_to_int_buffer);
                    Base::add_to_output(&space, 1);
                    Base::add_to_output(string_to_int_buffer, len);
                }
            }
            Base::add_to_output(&newline, 1);

            // Reset counts
            for(int64_t idx : nonzero_count_indices){
                counts[idx] = 0;
            }
            nonzero_count_indices.clear();
        }
    }

    virtual void finish(){
        Base::out->write(Base::output_buffer.data(), Base::output_buffer.size());
    }

    virtual ~ThresholdPseudoaligner() {} 
};

template<class coloring_t>
class IntersectionPseudoaligner : public DispatcherConsumerCallback, Pseudoaligner_Base<coloring_t>{

    private:

    typedef Pseudoaligner_Base<coloring_t> Base;

    IntersectionPseudoaligner(const IntersectionPseudoaligner&); // No copying
    IntersectionPseudoaligner& operator=(const IntersectionPseudoaligner&); // No copying

    public:

        IntersectionPseudoaligner(const plain_matrix_sbwt_t* SBWT, const coloring_t* coloring, ParallelBaseWriter* out, bool reverse_complements, int64_t output_buffer_capacity) : Pseudoaligner_Base<coloring_t>(SBWT, coloring, out, reverse_complements, output_buffer_capacity){}

        // Returns the color set
        vector<int64_t> do_intersections_on_color_id_buffers_with_reverse_complements(){
            int64_t n_kmers = Base::color_set_id_buffer.size();

            bool first_nonempty_union_found = false;
            typename coloring_t::colorset_type result;
            for(int64_t i = 0; i < n_kmers; i++){
                if(i > 0
                && (Base::color_set_id_buffer[i] == Base::color_set_id_buffer[i-1])
                && (Base::rc_color_set_id_buffer[n_kmers-1-i] == Base::rc_color_set_id_buffer[n_kmers-1-i+1])){
                    continue; // This pair of color set ids was already intersected in the previous iteration
                }

                int64_t fw_id = Base::color_set_id_buffer[i];
                int64_t rc_id = Base::rc_color_set_id_buffer[n_kmers-1-i];

                typename coloring_t::colorset_type cs; // For union. TODO: Reuse a union buffer
                if(fw_id == -1 && rc_id == -1) continue; // Neither direction is found
                else if(fw_id == -1 && rc_id >= 0) cs = Base::coloring->get_color_set_by_color_set_id(rc_id);
                else if(fw_id >= 0 && rc_id == -1) cs = Base::coloring->get_color_set_by_color_set_id(fw_id);
                else if(fw_id >= 0 && rc_id >= 0){
                    // Take union of forward and reverse complement
                    cs = Base::coloring->get_color_set_by_color_set_id(fw_id);
                    cs.do_union(Base::coloring->get_color_set_by_color_set_id(rc_id));
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
            int64_t n_kmers = Base::color_set_id_buffer.size();

            bool first_nonempty_color_set_found = false;
            typename coloring_t::colorset_type result;
            for(int64_t i = 0; i < n_kmers; i++){
                if(i > 0  && (Base::color_set_id_buffer[i] == Base::color_set_id_buffer[i-1])){
                    continue; // This color set was already intersected in the previous iteration
                }
                if(Base::color_set_id_buffer[i] == -1) continue; // k-mer not found

                const typename coloring_t::colorset_type& cs = Base::coloring->get_color_set_by_color_set_id(Base::color_set_id_buffer[i]);
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

        virtual void callback(const char* S, int64_t S_size, int64_t string_id, std::array<uint8_t, 8> metadata){

            (void) metadata; // Ignored

            char string_to_int_buffer[32]; // Enough space for a 64-bit integer in ascii
            char newline = '\n';
            char space = ' ';

            Base::color_set_id_buffer.resize(0);
            Base::rc_color_set_id_buffer.resize(0);
            // Clearing the buffers like this might look bad for performance at first glance because
            // we will then need to allocate new space for new elements that will be pushed to the buffers.
            // But in fact it's ok because resize is not supposed to affect the internal capacity of the vector.
            // cppreference.com says:
            //     "Vector capacity is never reduced when resizing to smaller size because that would
            //      invalidate all iterators, rather than only the ones that would be invalidated by the
            //      equivalent sequence of pop_back() calls."

            if(S_size < Base::k){
                write_log("Warning: query is shorter than k", LogLevel::MINOR);
                int64_t len = fast_int_to_string(string_id, string_to_int_buffer);
                Base::add_to_output(string_to_int_buffer, len);
                Base::add_to_output(&newline, 1);
            }
            else{
                vector<int64_t> colex_ranks = Base::SBWT->streaming_search(S, S_size); // TODO: version that pushes to existing buffer?
                Base::push_color_set_ids_to_buffer(colex_ranks, Base::color_set_id_buffer);
                if(Base::reverse_complements){
                    Base::push_color_set_ids_to_buffer(Base::get_rc_colex_ranks(S, S_size), Base::rc_color_set_id_buffer);
                }

                vector<int64_t> intersection;
                if(Base::reverse_complements) intersection = do_intersections_on_color_id_buffers_with_reverse_complements();
                else intersection = do_intersections_on_color_id_buffers_without_reverse_complements();

                int64_t len = fast_int_to_string(string_id, string_to_int_buffer);
                Base::add_to_output(string_to_int_buffer, len);
                for(color_t x : intersection){
                    len = fast_int_to_string(x, string_to_int_buffer);
                    Base::add_to_output(&space, 1);
                    Base::add_to_output(string_to_int_buffer, len);
                }
                Base::add_to_output(&newline, 1);
            }
        }

        virtual void finish(){
            Base::out->write(Base::output_buffer.data(), Base::output_buffer.size());
        }
};

template<typename instream_t, typename outstream_t>
void sort_parallel_output_file(instream_t& instream, outstream_t& outstream){
    set<pair<int64_t, string> > Q; // Priority queue with pairs (priority, content)
    string line;
    vector<string> tokens;
    int64_t current_query_id = 0;

    while(getline(instream,line)){
        stringstream ss(line);
        int64_t priority; ss >> priority;
        Q.insert({priority, line + "\n"});
        while(Q.begin()->first == current_query_id){
            outstream << Q.begin()->second;
            Q.erase(Q.begin());
            current_query_id++;
        }
    }
    assert(Q.empty());
}

// If outfile is empty, creates a writer to cout
std::unique_ptr<ParallelBaseWriter> create_writer(const string& outfile, bool gzipped);

void call_sort_parallel_output_file(const string& outfile, bool gzipped);

// Metadata stream for non-existent metadata
class Null_Metadata_Stream : public Metadata_Stream{

private:

std::array<uint8_t, 8> dummy;

public:
    virtual std::array<uint8_t, 8> next(){
        return dummy;
    }
};

template<typename coloring_t, typename sequence_reader_t>
void pseudoalign_intersected(const plain_matrix_sbwt_t& SBWT, const coloring_t& coloring, int64_t n_threads, sequence_reader_t& reader, std::string outfile, bool reverse_complements, int64_t buffer_size, bool gzipped, bool sorted_output){

    std::unique_ptr<ParallelBaseWriter> out = create_writer(outfile, gzipped);

    vector<DispatcherConsumerCallback*> threads;
    for (int64_t i = 0; i < n_threads; i++) {
        IntersectionPseudoaligner<coloring_t>* T = new IntersectionPseudoaligner<coloring_t>(&SBWT, &coloring, out.get(), reverse_complements, buffer_size);
        threads.push_back(T);
    }

    Null_Metadata_Stream metadata_stream; // No metadata
    run_dispatcher(threads, reader, &metadata_stream, buffer_size);

    // Clean up
    for (DispatcherConsumerCallback* t : threads) delete t;
    out->flush();
    out.reset(); // Free the memory

    if (sorted_output) call_sort_parallel_output_file(outfile, gzipped);
    
}

template<typename coloring_t, typename sequence_reader_t>
void pseudoalign_thresholded(const plain_matrix_sbwt_t& SBWT, const coloring_t& coloring, int64_t n_threads, sequence_reader_t& reader, std::string outfile, bool reverse_complements, int64_t buffer_size, bool gzipped, bool sorted_output, double threshold, bool ignore_unknown){

    std::unique_ptr<ParallelBaseWriter> out = create_writer(outfile, gzipped);

    vector<DispatcherConsumerCallback*> threads;
    for (int64_t i = 0; i < n_threads; i++) {
        ThresholdPseudoaligner<coloring_t>* T = new ThresholdPseudoaligner<coloring_t>(&SBWT, &coloring, out.get(), reverse_complements, buffer_size, threshold, ignore_unknown);
        threads.push_back(T);
    }

    Null_Metadata_Stream metadata_stream; // No metadata
    run_dispatcher(threads, reader, &metadata_stream, buffer_size);

    // Clean up
    for (DispatcherConsumerCallback* t : threads) delete t;
    out->flush();
    out.reset(); // Free the memory

    if (sorted_output) call_sort_parallel_output_file(outfile, gzipped);
    
}