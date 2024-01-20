#pragma once

#include <string>
#include "sbwt/SBWT.hh"
#include "coloring/Coloring.hh"
#include "SeqIO/SeqIO.hh"
#include "ThreadPool.hh"
#include "variants.hh"

using namespace std;
using namespace sbwt;

typedef int64_t color_t;

// returns a vector where element i is the ref ids aligned with query i
vector<vector<int64_t> > parse_pseudoalignment_output_format_from_disk(string filename);

// If outfile is empty, creates a writer to cout
std::unique_ptr<ParallelBaseWriter> create_writer(const string& outfile, bool gzipped);

void call_sort_parallel_output_file(const string& outfile, bool gzipped);

namespace pseudoalignment{ // Helper classes for pseudoalignment.

template<class coloring_t>
class Pseudoaligner_Base{

    // This class contains code that is common to both the intersection and threshold pseudualigners

private:

    // Flushes at the next newline when output_buffer_flush_threshold is exceeded
    template<typename output_stream_t>
    void write_out(const char* data, int64_t data_length, output_stream_t& output_writer, vector<char>& buffer){
        for(int64_t i = 0; i < data_length; i++){
            buffer.push_back(data[i]);
            if(buffer.size() > output_buffer_flush_threshold && data[i] == '\n'){
                output_writer->write(buffer.data(), buffer.size());
                *total_bytes_written += buffer.size();
                buffer.clear(); // Let's hope this keeps the reserved capacity of the vector intact
            }
        }
    }

public:

    const plain_matrix_sbwt_t* SBWT; // Not owned by this class
    const coloring_t* coloring; // Not owned by this class
    shared_ptr<ParallelBaseWriter> out;
    optional<shared_ptr<ParallelBaseWriter>> aux_out;

    bool reverse_complements;
    int64_t k;
    double relevant_kmers_fraction;
    bool sort_hits;

    // Buffer for reverse-complementing strings
    vector<char> rc_buffer;

    // Buffers for printing. We want to have local buffers for each thread to avoid having to call the
    // parallel writer so often to avoid locking the writers from other threads.
    vector<char> output_buffer;
    vector<char> aux_output_buffer;
    int64_t output_buffer_flush_threshold; // Used also as a threshold for aux output buffer

    // Buffer for storing color set ids
    vector<int64_t> color_set_id_buffer;
    vector<int64_t> rc_color_set_id_buffer;

    // Statistics to print. These will be read by a printer thread while
    // they are modified, so they need to be atomic.
    atomic<int64_t>* total_length_of_sequence_processed;
    atomic<int64_t>* total_bytes_written;

    // For output writing
    char int_to_string_buffer[32]; // Enough space for a 64-bit integer in ascii

    Pseudoaligner_Base(const plain_matrix_sbwt_t* SBWT, const coloring_t* coloring, shared_ptr<ParallelBaseWriter> out, optional<shared_ptr<ParallelBaseWriter>> aux_out, bool reverse_complements, int64_t output_buffer_capacity, atomic<int64_t>* total_length_of_sequence_processed, atomic<int64_t>* total_bytes_written, double relevant_kmers_fraction, bool sort_hits){
        this->SBWT = SBWT;
        this->coloring = coloring;
        this->out = out;
        this->aux_out = aux_out;
        this->reverse_complements = reverse_complements;
        this->output_buffer_flush_threshold = output_buffer_capacity;
        this->k = SBWT->get_k();
        this->total_length_of_sequence_processed = total_length_of_sequence_processed;
        this->total_bytes_written = total_bytes_written;
        this->relevant_kmers_fraction = relevant_kmers_fraction;
        this->sort_hits = sort_hits;
        rc_buffer.resize(1 << 10); // 1 kb. Will be resized if needed
        output_buffer.reserve(output_buffer_capacity);
        aux_output_buffer.reserve(output_buffer_capacity);
    }

    // Flushes at the next newline when output_buffer_flush_threshold is exceeded
    void add_to_output(const char* data, int64_t data_length){
        write_out(data, data_length, out, output_buffer);
    }

    void add_to_aux_output(const char* data, int64_t data_length){
        if(aux_out.has_value()) write_out(data, data_length, aux_out.value(), aux_output_buffer);
    }

    // Hit counts is such that hit_counts[color] is the number of hits of color
    void report_results_for_seq(const char* seq, int64_t seq_id, int64_t seq_len, vector<int64_t>& hit_colors, const optional<vector<int64_t>>& hit_counts, int64_t n_kmers_found_in_index){
        if(sort_hits) std::sort(hit_colors.begin(), hit_colors.end()); 

        // Add the sequence id
        int64_t len = fast_int_to_string(seq_id, int_to_string_buffer);
        add_to_output(int_to_string_buffer, len);

        // Add the color hits
        for(color_t x : hit_colors){
            len = fast_int_to_string(x, int_to_string_buffer);
            add_to_output(" ", 1);
            add_to_output(int_to_string_buffer, len);
        }
        add_to_output("\n", 1);

        // Write the aux info, if needed
        if(this->aux_out.has_value()){
            // Sequence id
            add_to_aux_output("{\"seq-id\":", 10);
            add_to_aux_output(int_to_string_buffer, fast_int_to_string(seq_id, int_to_string_buffer));

            // Number of kmers
            int64_t total_kmers = max((int64_t)0, seq_len - k + 1);
            add_to_aux_output(",\"kmers\":", 9);
            add_to_aux_output(int_to_string_buffer, fast_int_to_string(total_kmers, int_to_string_buffer));

            // Number of relevant k-mers
            add_to_aux_output(",\"relevant\":", 12);
            add_to_aux_output(int_to_string_buffer, fast_int_to_string(n_kmers_found_in_index, int_to_string_buffer));

            // Matching statistics
            add_to_aux_output(",\"MS\":[", 7);
            for(int64_t i = 0; i < seq_len - k + 1; i++){
                pair<pair<int64_t, int64_t>, int64_t> result = this->SBWT->partial_search(seq + i, k);
                int64_t match_length = result.second;
                if(i > 0) add_to_aux_output(",", 1);
                add_to_aux_output(int_to_string_buffer, fast_int_to_string(match_length, int_to_string_buffer));
            }
            add_to_aux_output("]", 1);

            // Hit counts
            if(hit_counts.has_value()){
                add_to_aux_output(",\"hit-counts\":{", 15);
                for(int64_t i = 0; i < hit_colors.size(); i++){
                    if(i > 0) add_to_aux_output(",", 1);
                    add_to_aux_output("\"", 1);
                    add_to_aux_output(int_to_string_buffer, fast_int_to_string(hit_colors[i], int_to_string_buffer));
                    add_to_aux_output("\":", 2);
                    add_to_aux_output(int_to_string_buffer, fast_int_to_string((*hit_counts)[hit_colors[i]], int_to_string_buffer));
                }
                add_to_aux_output("}", 1);
            }

            add_to_aux_output("}\n", 2);
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

    ~Pseudoaligner_Base(){
        // Flush remaining output
        out->write(output_buffer.data(), output_buffer.size());
        *total_bytes_written += output_buffer.size();

        if(aux_out.has_value()){
            aux_out.value()->write(aux_output_buffer.data(), aux_output_buffer.size());
            *total_bytes_written += aux_output_buffer.size();
        }
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

class WorkBatch{
    public:
        unique_ptr<vector<char>> seqs_concat;
        unique_ptr<vector<int64_t>> starts; // Has an end sentinel one past the last one
        unique_ptr<vector<int64_t>> seq_ids;

    // Default constructor
    WorkBatch(){
        seqs_concat = make_unique<vector<char>>();
        starts = make_unique<vector<int64_t>>();
        seq_ids = make_unique<vector<int64_t>>();
    }

    // Move constructor
    WorkBatch(WorkBatch&& other){

        // Move the stuff in
        this->seqs_concat = move(other.seqs_concat);
        this->starts = move(other.starts);
        this->seq_ids = move(other.seq_ids);

        // Clear the other
        other.seqs_concat = make_unique<vector<char>>();
        other.starts = make_unique<vector<int64_t>>();
        other.seq_ids = make_unique<vector<int64_t>>();
    }
};


// Context for the worker threads
template<typename coloring_t>
struct WorkerContext{
    const plain_matrix_sbwt_t* SBWT;
    const coloring_t* coloring;
    const bool reverse_complements;
    const double threshold;
    const bool ignore_unknown;
    const bool sort_hits;
    const int64_t output_buffer_size;
    
    // Statistics to print. These will be read by a printer thread while
    // they are modified, so they need to be atomic.
    atomic<int64_t>* total_length_of_sequence_processed = 0;
    atomic<int64_t>* total_bytes_written = 0;

    shared_ptr<ParallelBaseWriter> writer;
    std::optional<shared_ptr<ParallelBaseWriter>> aux_writer;

    // Don't move these fields up because we're using the struct initializer list syntax which depends on the order
    double relevant_kmers_fraction;

};

template <typename coloring_t>
class ThresholdWorker : public BaseWorkerThread<WorkBatch>, Pseudoaligner_Base<coloring_t>{
    public:

    typedef WorkBatch work_item_t;
    typedef WorkerContext<coloring_t> Context;
    typedef Pseudoaligner_Base<coloring_t> Base;

    double count_threshold; // Fraction of k-mers that need to be found to report pseudoalignment to a color
    bool ignore_unknown_kmers = false; // Ignore k-mers that do not exist in the de Bruijn graph or have no colors


    // State used during callback
    vector<int64_t> counts; // counts[i] = number of occurrences of color i
    vector<int64_t> nonzero_count_indices; // Indices in this->counts that have a non-zero value
    vector<int64_t> color_buffer; // Reused buffer for storing colors
    vector<int64_t> hits; // Pseudoalignment hits to report

    ThresholdWorker(WorkerContext<coloring_t> context) :
        Pseudoaligner_Base<coloring_t>(context.SBWT, context.coloring, context.writer, context.aux_writer, context.reverse_complements, context.output_buffer_size, context.total_length_of_sequence_processed, context.total_bytes_written, context.relevant_kmers_fraction, context.sort_hits), count_threshold(context.threshold), ignore_unknown_kmers(context.ignore_unknown){
        counts.resize(context.coloring->largest_color() + 1); // Initializes counts to zeroes
    }


    void process_sequence(const char* S, int64_t S_size, int64_t string_id){

        Base::color_set_id_buffer.resize(0);
        Base::rc_color_set_id_buffer.resize(0);

        if(S_size < Base::k){
            write_log("Warning: query is shorter than k", LogLevel::MINOR);
            hits.clear();
            Base::report_results_for_seq(S, string_id, S_size, hits, counts, 0);
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
                    color_buffer.clear();
                    fw_set.push_colors_to_vector(color_buffer);
                    for(int64_t color : color_buffer){
                        has_at_least_one_color = true;
                        if(counts[color] == 0) nonzero_count_indices.push_back(color);
                        counts[color] += run_length;
                    }

                    n_kmers_with_at_least_1_color += has_at_least_one_color * run_length;

                    run_length = 0; // Reset the run
                }
            }

            // Print the colors of all counters that are above threshold and the relevant k-mers fraction
            hits.clear();
            for(int64_t color : nonzero_count_indices){
                int64_t count = counts[color];
                int64_t effective_kmers = ignore_unknown_kmers ? n_kmers_with_at_least_1_color : n_kmers;
                if(count >= effective_kmers * count_threshold && (double)effective_kmers / n_kmers >= Base::relevant_kmers_fraction){
                    // Add to list of reported colors
                    hits.push_back(color);
                }
            }
            Base::report_results_for_seq(S, string_id, S_size, hits, counts, n_kmers_with_at_least_1_color);

            // Reset counts
            for(int64_t idx : nonzero_count_indices){
                counts[idx] = 0;
            }
            nonzero_count_indices.clear();
        }
    }

    // This function should only use local variables and protected shared variables
    virtual void process_work_item(WorkBatch item){
        for(int64_t i = 0; i < (int64_t)item.starts->size() - 1; i++){
            int64_t start = (*item.starts)[i];
            int64_t end = (*item.starts)[i+1];
            int64_t seq_id = (*item.seq_ids)[i];
            process_sequence(item.seqs_concat->data() + start, end-start, seq_id);
            *Base::total_length_of_sequence_processed += end - start;
        }
    }

    // This function is called after every work item. It is called so
    // that only one thread at a time is executing the function.
    virtual void critical_section(){
        // We have no critical section because the output synchronization is handled by the output writer
    }
};

template <typename coloring_t>
class IntersectionWorker : public BaseWorkerThread<WorkBatch>, Pseudoaligner_Base<coloring_t>{
    public:

    typedef WorkBatch work_item_t;
    typedef WorkerContext<coloring_t> Context;
    typedef Pseudoaligner_Base<coloring_t> Base;

    IntersectionWorker(WorkerContext<coloring_t> context) :
        Pseudoaligner_Base<coloring_t>(context.SBWT, context.coloring, context.writer, context.aux_writer, context.reverse_complements, context.output_buffer_size, context.total_length_of_sequence_processed, context.total_bytes_written, context.relevant_kmers_fraction, context.sort_hits){}

    void process_sequence(const char* S, int64_t S_size, int64_t string_id){

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
            vector<color_t> empty;
            Base::report_results_for_seq(S, string_id, S_size, empty, nullopt, 0); // Intersection does not compute hit counts, so we pass a nullopt
        }
        else{
            vector<int64_t> colex_ranks = Base::SBWT->streaming_search(S, S_size); // TODO: version that pushes to existing buffer?
            Base::push_color_set_ids_to_buffer(colex_ranks, Base::color_set_id_buffer);
            if(Base::reverse_complements){
                Base::push_color_set_ids_to_buffer(Base::get_rc_colex_ranks(S, S_size), Base::rc_color_set_id_buffer);
            }

            vector<int64_t> intersection; int64_t n_nonempty;
            if(Base::reverse_complements) tie(intersection, n_nonempty) = do_intersections_on_color_id_buffers_with_reverse_complements();
            else tie(intersection, n_nonempty) = do_intersections_on_color_id_buffers_without_reverse_complements();

            if((double)n_nonempty / (S_size - Base::k + 1) >= Base::relevant_kmers_fraction)
                Base::report_results_for_seq(S, string_id, S_size, intersection, nullopt, n_nonempty); // Intersection does not compute hit counts, so we pass a nullopt
        }
    }

    // Returns the color set and the number of non-empty colorsets in the query
    pair<vector<int64_t>, int64_t> do_intersections_on_color_id_buffers_with_reverse_complements(){
        int64_t n_kmers = Base::color_set_id_buffer.size();

        int64_t prev_colorset_size = 0;
        int64_t n_nonempty = 0;
        
        typename coloring_t::colorset_type result;
        for(int64_t i = 0; i < n_kmers; i++){
            if(i > 0
            && (Base::color_set_id_buffer[i] == Base::color_set_id_buffer[i-1])
            && (Base::rc_color_set_id_buffer[n_kmers-1-i] == Base::rc_color_set_id_buffer[n_kmers-1-i+1])){
                if(prev_colorset_size > 0) n_nonempty++;
                continue; // This pair of color set ids was already intersected in the previous iteration
            }

            int64_t fw_id = Base::color_set_id_buffer[i];
            int64_t rc_id = Base::rc_color_set_id_buffer[n_kmers-1-i];

            // Figure out the color set
            typename coloring_t::colorset_type cs; // For union. TODO: Reuse a union buffer
            if(fw_id == -1 && rc_id == -1){
                prev_colorset_size = 0;
                continue; // Neither direction is found
            }
            else if(fw_id == -1 && rc_id >= 0) cs = Base::coloring->get_color_set_by_color_set_id(rc_id);
            else if(fw_id >= 0 && rc_id == -1) cs = Base::coloring->get_color_set_by_color_set_id(fw_id);
            else if(fw_id >= 0 && rc_id >= 0){
                // Take union of forward and reverse complement
                cs = Base::coloring->get_color_set_by_color_set_id(fw_id);
                cs.do_union(Base::coloring->get_color_set_by_color_set_id(rc_id));
            }

            if(cs.size() > 0){
                if(n_nonempty == 0) result = cs; // This is the first nonempty color set
                else result.intersection(cs); // Intersection
                n_nonempty++;
            }
            prev_colorset_size = cs.size();
        }
        return {result.get_colors_as_vector(), n_nonempty};
    }

    // Returns the color set and the number of non-empty colorsets in the query
    pair<vector<int64_t>, int64_t> do_intersections_on_color_id_buffers_without_reverse_complements(){
        int64_t n_kmers = Base::color_set_id_buffer.size();

        int64_t prev_colorset_size = 0;
        int64_t n_nonempty = 0;
        typename coloring_t::colorset_type result;
        for(int64_t i = 0; i < n_kmers; i++){
            if(i > 0  && (Base::color_set_id_buffer[i] == Base::color_set_id_buffer[i-1])){
                // This color set was already intersected in the previous iteration
                if(prev_colorset_size > 0) n_nonempty++;
                // prev_colorset_size unchanged
            } else if(Base::color_set_id_buffer[i] == -1){
                // k-mer not found
                prev_colorset_size = 0;
            } else{
                // k-mer is found and it has a different color set from the previous one
                const typename coloring_t::colorset_type::view_t cs = Base::coloring->get_color_set_by_color_set_id(Base::color_set_id_buffer[i]);
                if(cs.size() > 0){
                    if(n_nonempty == 0) result = cs; // This is the first nonempty color set
                    else result.intersection(cs); // Intersection
                    n_nonempty++;
                }
                prev_colorset_size = cs.size();
            }
        }
        return {result.get_colors_as_vector(), n_nonempty};
    }

    // This function should only use local variables and protected shared variables
    virtual void process_work_item(WorkBatch item){
        for(int64_t i = 0; i < (int64_t)item.starts->size() - 1; i++){
            int64_t start = (*item.starts)[i];
            int64_t end = (*item.starts)[i+1];
            int64_t seq_id = (*item.seq_ids)[i];
            process_sequence(item.seqs_concat->data() + start, end-start, seq_id);
            *Base::total_length_of_sequence_processed += end - start;
        }
    }

    // This function is called after every work item. It is called so
    // that only one thread at a time is executing the function.
    virtual void critical_section(){
        // We have no critical section because the output synchronization is handled by the output writer
    }
};

// Forwards work to either the an intersection worker or a threshold worker
template <typename coloring_t>
class Worker : public BaseWorkerThread<WorkBatch>{
    public:

        unique_ptr<BaseWorkerThread<WorkBatch>> inner_worker;

        Worker(WorkerContext<coloring_t> context){
            // Initialize the correct inner worker
            // If there is an aux writer then we can't use the intersection worker because
            // It does not keep track of the counts which should be reported in the aux file.
            if(context.threshold == 1 && !(context.aux_writer.has_value()))
                inner_worker = make_unique<IntersectionWorker<coloring_t>>(context);
            else
                inner_worker = make_unique<ThresholdWorker<coloring_t>>(context);
        }

        // This function should only use local variables and protected shared variables
        virtual void process_work_item(WorkBatch item){
            inner_worker->process_work_item(move(item));
        }

        // This function is called after every work item. It is called so
        // that only one thread at a time is executing the function.
        virtual void critical_section(){
            inner_worker->critical_section();
        }

        virtual ~Worker() = default;
};

void print_thread(atomic<int64_t>* total_length_of_sequence_processed, atomic<int64_t>* total_bytes_written, atomic<bool>* stop_printing);

template<typename sequence_reader_t, typename coloring_t>
void push_work_batches(int64_t buffer_size, sequence_reader_t& reader, ThreadPool<Worker<coloring_t>, pseudoalignment::WorkBatch>& TP){
    // Start creating work batches
    int64_t batch_push_threshold = buffer_size; // A batch is pushed to the thread pool when it reaches this size
    WorkBatch wb;

    int64_t seq_id = 0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;

        // Add the read to the batch
        wb.starts->push_back(wb.seqs_concat->size());
        wb.seq_ids->push_back(seq_id);
        for(int64_t i = 0 ; i < len; i++){
            wb.seqs_concat->push_back(reader.read_buf[i]);
        }

        if(wb.seqs_concat->size() >= batch_push_threshold){
            // Push the batch
            wb.starts->push_back(wb.seqs_concat->size()); // End sentinel
            TP.add_work(std::move(wb), wb.seqs_concat->size());
            // Moving the batch also clears it
        }

        seq_id++;
    }

    // Push the last batch
    if(wb.seqs_concat->size() > 0){
        wb.starts->push_back(wb.seqs_concat->size()); // End sentinel
        TP.add_work(std::move(wb), wb.seqs_concat->size());
    }

}

} // End namespace pseudoalignment

template<typename coloring_t, typename sequence_reader_t>
void pseudoalign(const plain_matrix_sbwt_t& SBWT, const coloring_t& coloring, int64_t n_threads, sequence_reader_t& reader, std::string outfile, std::string aux_info_file, bool reverse_complements, int64_t buffer_size, bool gzipped, bool sorted_output, double threshold, bool ignore_unknown, double relevant_kmers_fraction, bool sort_hits){

    using namespace pseudoalignment;

    // Set up context (= commmon variables for all workers).

    std::shared_ptr<ParallelBaseWriter> out = create_writer(outfile, gzipped);

    std::optional<std::shared_ptr<ParallelBaseWriter>> aux_out;
    if(aux_info_file != "") aux_out = create_writer(aux_info_file, gzipped);

    atomic<int64_t> total_length_of_sequence_processed = 0; // For printing progress
    atomic<int64_t> total_bytes_written = 0; // For printing progress
    WorkerContext<coloring_t> context = {&SBWT, &coloring, reverse_complements, threshold, ignore_unknown, sort_hits, buffer_size, &total_length_of_sequence_processed, &total_bytes_written, out, aux_out, relevant_kmers_fraction};

    // Create workers
    vector<unique_ptr<Worker<coloring_t>>> workers;
    vector<Worker<coloring_t>*> worker_ptrs;
    for(int64_t i = 0; i < n_threads; i++){
        workers.push_back(make_unique<Worker<coloring_t>>(context));
        worker_ptrs.push_back(workers.back().get());
    }

    // Launch a thread that prints progress every second until done
    atomic<bool> stop_printing = false; // The thread will stop when this is set to true
    std::thread print_thread(pseudoalignment::print_thread, &total_length_of_sequence_processed, &total_bytes_written, &stop_printing);

    // Create a worker thread pool
    ThreadPool<Worker<coloring_t>, WorkBatch> TP(worker_ptrs, buffer_size);

    try{ // For some reason exceptions are not propagates up to main from here, so we catch them and terminate the program here
        push_work_batches(buffer_size, reader, TP);
    } catch (const std::runtime_error &e){
        std::cerr << "Runtime error: " << e.what() << '\n';
        std::terminate();
    }

    TP.join_threads();
    workers.clear(); // This will delete the workers, which will flush their internal buffers to the common output buffer

    // Flush the common output buffers
    out->flush();
    if(aux_out.has_value()) aux_out.value()->flush();

    // Terminate the print thread
    stop_printing = true;
    print_thread.join();

    if (sorted_output) call_sort_parallel_output_file(outfile, gzipped);
    if (sorted_output && aux_out.has_value()) call_sort_parallel_output_file(aux_info_file, gzipped);
    
}