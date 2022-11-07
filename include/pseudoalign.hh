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

        Pseudoaligner_Base(const plain_matrix_sbwt_t* SBWT, const coloring_t* coloring, ParallelBaseWriter* out, bool reverse_complements, LL output_buffer_capacity){
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

};

template<class coloring_t>
class ThresholdPseudoaligner : public DispatcherConsumerCallback, Pseudoaligner_Base<coloring_t>{

private:

typedef Pseudoaligner_Base<coloring_t> Base;
double count_threshold; // Fraction of k-mers that need to be found to report pseudoalignment to a color

public:

    ThresholdPseudoaligner(const plain_matrix_sbwt_t* SBWT, const coloring_t* coloring, ParallelBaseWriter* out, bool reverse_complements, LL output_buffer_capacity, double count_threshold) : Pseudoaligner_Base<coloring_t>(SBWT, coloring, out, reverse_complements, output_buffer_capacity), count_threshold(count_threshold){}

    virtual void callback(const char* S, LL S_size, int64_t string_id){
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
        }
        else{
            vector<int64_t> colex_ranks = Base::SBWT->streaming_search(S, S_size); // TODO: version that pushes to existing buffer?
            Base::push_color_set_ids_to_buffer(colex_ranks, Base::color_set_id_buffer);
            vector<int64_t> rc_colex_ranks;
            if(Base::reverse_complements){
                while(S_size > Base::rc_buffer.size()){
                    Base::rc_buffer.resize(Base::rc_buffer.size()*2);
                }
                memcpy(Base::rc_buffer.data(), S, S_size);
                reverse_complement_c_string(Base::rc_buffer.data(), S_size); // There is no null at the end but that is ok
                rc_colex_ranks = Base::SBWT->streaming_search(Base::rc_buffer.data(), S_size);
                Base::push_color_set_ids_to_buffer(rc_colex_ranks, Base::rc_color_set_id_buffer);
            }

            // Todo: use something better than a std::map for the counters
            unordered_map<int64_t, int64_t> counts; // color id -> count of that color id
            typename coloring_t::colorset_type fw_set;
            int64_t n_kmers = S_size - Base::k  + 1;
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
                        typename coloring_t::colorset_type::view_t rc_set_view = Base::coloring->get_color_set_by_color_set_id(Base::rc_color_set_id_buffer[S_size - Base::k  - kmer_idx]);
                        fw_set.do_union(rc_set_view);
                    }

                    // Add the run length to the counts
                    for(int64_t color : fw_set.get_colors_as_vector()){
                        counts[color] += run_length;
                    }

                    run_length = 0; // Reset the run
                }
            }

            // Print the colors of all counters that are above threshold
            int64_t len = fast_int_to_string(string_id, string_to_int_buffer);
            Base::add_to_output(string_to_int_buffer, len); // String id
            for(auto [color, count] : counts){
                if(count >= (S_size - Base::k + 1) * count_threshold){
                    // Report color
                    len = fast_int_to_string(color, string_to_int_buffer);
                    Base::add_to_output(&space, 1);
                    Base::add_to_output(string_to_int_buffer, len);
                }
            }
            Base::add_to_output(&newline, 1);
        }
    }

    virtual void finish(){
        Base::out->write(Base::output_buffer.data(), Base::output_buffer_size);
        Base::output_buffer_size = 0;
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

        IntersectionPseudoaligner(const plain_matrix_sbwt_t* SBWT, const coloring_t* coloring, ParallelBaseWriter* out, bool reverse_complements, LL output_buffer_capacity) : Pseudoaligner_Base<coloring_t>(SBWT, coloring, out, reverse_complements, output_buffer_capacity){}

        // Returns the color set
        vector<int64_t> do_intersections_on_color_id_buffers_with_reverse_complements(){
            LL n_kmers = Base::color_set_id_buffer.size();

            bool first_nonempty_union_found = false;
            typename coloring_t::colorset_type result;
            for(LL i = 0; i < n_kmers; i++){
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
            LL n_kmers = Base::color_set_id_buffer.size();

            bool first_nonempty_color_set_found = false;
            typename coloring_t::colorset_type result;
            for(LL i = 0; i < n_kmers; i++){
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

        virtual void callback(const char* S, LL S_size, int64_t string_id){
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
                vector<int64_t> rc_colex_ranks;
                if(Base::reverse_complements){
                    while(S_size > Base::rc_buffer.size()){
                        Base::rc_buffer.resize(Base::rc_buffer.size()*2);
                    }
                    memcpy(Base::rc_buffer.data(), S, S_size);
                    reverse_complement_c_string(Base::rc_buffer.data(), S_size); // There is no null at the end but that is ok
                    rc_colex_ranks = Base::SBWT->streaming_search(Base::rc_buffer.data(), S_size);
                    Base::push_color_set_ids_to_buffer(rc_colex_ranks, Base::rc_color_set_id_buffer);
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
            Base::out->write(Base::output_buffer.data(), Base::output_buffer_size);
            Base::output_buffer_size = 0;
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

// If outfile is empty, creates a writer to cout
std::unique_ptr<ParallelBaseWriter> create_writer(const string& outfile, bool gzipped);

void call_sort_parallel_output_file(const string& outfile, bool gzipped);

template<typename coloring_t, typename sequence_reader_t>
void pseudoalign_intersected(const plain_matrix_sbwt_t& SBWT, const coloring_t& coloring, int64_t n_threads, sequence_reader_t& reader, std::string outfile, bool reverse_complements, int64_t buffer_size, bool gzipped, bool sorted_output){

    std::unique_ptr<ParallelBaseWriter> out = create_writer(outfile, gzipped);

    vector<DispatcherConsumerCallback*> threads;
    for (LL i = 0; i < n_threads; i++) {
        IntersectionPseudoaligner<coloring_t>* T = new IntersectionPseudoaligner<coloring_t>(&SBWT, &coloring, out.get(), reverse_complements, buffer_size);
        threads.push_back(T);
    }

    run_dispatcher(threads, reader, buffer_size);

    // Clean up
    for (DispatcherConsumerCallback* t : threads) delete t;
    out->flush();
    out.reset(); // Free the memory

    if (sorted_output) call_sort_parallel_output_file(outfile, gzipped);
    
}

template<typename coloring_t, typename sequence_reader_t>
void pseudoalign_thresholded(const plain_matrix_sbwt_t& SBWT, const coloring_t& coloring, int64_t n_threads, sequence_reader_t& reader, std::string outfile, bool reverse_complements, int64_t buffer_size, bool gzipped, bool sorted_output, double threshold){

    std::unique_ptr<ParallelBaseWriter> out = create_writer(outfile, gzipped);

    vector<DispatcherConsumerCallback*> threads;
    for (LL i = 0; i < n_threads; i++) {
        ThresholdPseudoaligner<coloring_t>* T = new ThresholdPseudoaligner<coloring_t>(&SBWT, &coloring, out.get(), reverse_complements, buffer_size, threshold);
        threads.push_back(T);
    }

    run_dispatcher(threads, reader, buffer_size);

    // Clean up
    for (DispatcherConsumerCallback* t : threads) delete t;
    out->flush();
    out.reset(); // Free the memory

    if (sorted_output) call_sort_parallel_output_file(outfile, gzipped);
    
}