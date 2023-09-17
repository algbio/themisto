#pragma once

#include <algorithm>
#include <string>
#include <utility>
#include <vector>
#include <cstdint>
#include <cstring>
#include <variant>
#include <mutex>

#include <sdsl/bit_vectors.hpp>

#include "core_kmer_marker.hh"
#include "backward_traversal.hh"
#include "Sparse_Uint_Array.hh"
#include "Coloring.hh"

#include "SeqIO/buffered_streams.hh"
#include "sbwt/EM_sort/bit_level_stuff.hh"
#include "sbwt/EM_sort/EM_sort.hh"
#include "sbwt/globals.hh"
#include "SeqIO/SeqIO.hh"

#include "WorkDispatcher.hh"
#include "hybrid_color_set.hh"
#include "Roaring_Color_Set.hh"
#include "Fixed_Width_Int_Color_Set.hh"

#include "ggcat.hh"

#include "ThreadPool.hh"
#include "Coloring_Builder.hh"

using namespace ggcat;

// An iterable GGCAT unitig database.
class GGCAT_unitig_database{

private:
    // No copying
    GGCAT_unitig_database(GGCAT_unitig_database const& other) = delete;
    GGCAT_unitig_database& operator=(GGCAT_unitig_database const& other) = delete;

public:

    GGCATInstance* instance;
    string graph_file;
    vector<string> color_names;
    int64_t k;

    GGCAT_unitig_database(vector<string>& filenames, int64_t mem_gigas, int64_t k, int64_t n_threads, bool canonical) : k(k) {

        GGCATConfig config;

        config.use_temp_dir = true;
        config.temp_dir = get_temp_file_manager().get_dir();
        config.memory = mem_gigas;
        config.prefer_memory = true;
        config.total_threads_count = n_threads;
        config.intermediate_compression_level = -1;

        config.use_stats_file = false;
        config.stats_file = "";

        // This leaks memory but it's only a few bytes. It can't be easily fixed because
        // this memory is allocated in the Rust API of GGCAT and freed only at the
        // end of the program.
        instance = GGCATInstance::create(config);
        
        graph_file = get_temp_file_manager().create_filename("",".fa");

        color_names.clear();
        for(int64_t i = 0; i < filenames.size(); i++){
            color_names.push_back(to_string(i));
        }

        std::string output_file = instance->build_graph_from_files(
            Slice<std::string>(filenames.data(), filenames.size()),
            graph_file,
            k,
            n_threads,
            !canonical,
            1,
            ExtraElaborationStep_UnitigLinks,
            true,
            Slice<std::string>(color_names.data(), color_names.size()),
            -1);

        vector<string> file_color_names = GGCATInstance::dump_colors(GGCATInstance::get_colormap_file(graph_file));

    }

    // The callback takes a unitig, the color set, and the is_same flag
    void iterate(std::function<void(const std::string&, const vector<int64_t>&, bool)> callback){

        vector<int64_t> prev_colors;
        std::mutex callback_mutex;

        auto outer_callback = [&](Slice<char> read, Slice<uint32_t> colors, bool same_colors){
            // Calls in callback provided by the caller of iterate.
            // WARNING: this function is called asynchronously from multiple threads, so it must be thread-safe.
            // Also the same_colors boolean is referred to the previous call of this function from the current thread.
            std::lock_guard<std::mutex> _lock(callback_mutex);
            try{
                string unitig = string(read.data, read.data + read.size);                
                if(same_colors){
                    callback(unitig, prev_colors, true);
                } else{
                    prev_colors.clear();
                    for (size_t i = 0; i < colors.size; i++){
                        prev_colors.push_back(colors.data[i]);
                    }
                    callback(unitig, prev_colors, false);
                }
            } catch(const std::exception& e){
                std::cerr << "Caught Error: " << e.what() << '\n';
                exit(1);
            }
        };

        this->instance->dump_unitigs(graph_file,k,1,true,outer_callback,true,-1);
    }

    string get_unitig_filename(){
        return graph_file;
    }

};

template<typename colorset_t = SDSL_Variant_Color_Set> 
requires Color_Set_Interface<colorset_t>
class Coloring_Builder_From_GGCAT{

    private:

    struct UnitigWorkBatch{
        vector<string> unitigs;
        vector<int64_t> color_set_ids;
    };

    class UnitigWorker : public BaseWorkerThread<UnitigWorkBatch>{

        public:

        struct context_t{

            // Read-only context
            const plain_matrix_sbwt_t* SBWT;
            const int64_t colorset_sampling_distance;

            // Mutable context
            Sparse_Uint_Array_Builder* builder; // Updates are applied here
            vector<pair<int64_t,int64_t>> color_set_pointer_updates; // Local buffer of updates
        };

        typedef UnitigWorkBatch work_item_t;

        private:

        context_t context;

        void queue_updates(const vector<int64_t>& colex_ranks, int64_t color_set_id){
            int64_t next_sample_counter = 0;
            for(int64_t i = (int64_t)colex_ranks.size()-1; i >= 0; i--){ // Iterate from end to start
                next_sample_counter++;

                if(colex_ranks[i] == -1){
                    next_sample_counter = 0;
                    continue;
                }

                if(i == (int64_t)colex_ranks.size() - 1 || next_sample_counter >= context.colorset_sampling_distance){
                    this->context.color_set_pointer_updates.push_back({colex_ranks[i], color_set_id});
                    next_sample_counter = 0;
                }
            }
        }

        public:


        UnitigWorker(context_t context) : context(context) {}


        // This function should only modify local variables and the updates buffer in the context
        virtual void process_work_item(UnitigWorkBatch batch){

            for(int64_t unitig_idx = 0; unitig_idx < batch.unitigs.size(); unitig_idx++){

                string& unitig = batch.unitigs[unitig_idx];
                int64_t color_set_id = batch.color_set_ids[unitig_idx];

                // Forward strand

                vector<int64_t> colex_ranks = context.SBWT->streaming_search(unitig);
                queue_updates(colex_ranks, color_set_id);

                // Reverse complement strand
                
                reverse_complement_c_string(&(unitig[0]), unitig.size());
                colex_ranks = context.SBWT->streaming_search(unitig);
                queue_updates(colex_ranks, color_set_id);

            }

        }

        // A mutex in the thread pool guards this function so that only one worker can be in this function at a time
        virtual void critical_section(){
            for(auto [node, pointer] : context.color_set_pointer_updates){
                context.builder->add(node, pointer);
            }
            context.color_set_pointer_updates.clear();
        }

    };

    void iterate_unitigs(GGCAT_unitig_database& unitig_database, Coloring<colorset_t>& coloring, ThreadPool<UnitigWorker, UnitigWorkBatch>& TP){
        // Create work batches for the workers
        UnitigWorkBatch work_batch;
        int64_t work_batch_max_size = 1000;
        int64_t color_set_id = -1;
        coloring.largest_color_id = -1;

        // This is the a callback given to the GGCAT unitig iterator.
        auto process_unitig_and_colors = [&](const string& unitig, const vector<int64_t>& colors, bool same_colors){
            // same_colors means that the color set of the current unitig is the same as the color set of
            // the previous.

            if(!same_colors){
                color_set_id++;
                // Keep track of maximum color
                for(int64_t x : colors) coloring.largest_color_id = max(x, coloring.largest_color_id);

                // Store color set
                coloring.sets.add_set(colors);
                coloring.total_color_set_length += colors.size();
            }

            work_batch.unitigs.push_back(unitig);
            work_batch.color_set_ids.push_back(color_set_id);

            if(work_batch.unitigs.size() == work_batch_max_size){
                TP.add_work(work_batch, 1);
                work_batch.unitigs.clear();
                work_batch.color_set_ids.clear();
            }
        };

        unitig_database.iterate(process_unitig_and_colors);

        if(work_batch.unitigs.size() > 0){ // Last batch
            TP.add_work(work_batch, 1);
        }

        TP.join_threads();
    }

    public:

    // Colored unitig stream database produce canonical bidirected unitigs
    void build_from_colored_unitigs(Coloring<colorset_t>& coloring,
                        const plain_matrix_sbwt_t& SBWT,
                        const std::int64_t ram_bytes,
                        const std::int64_t n_threads,
                        int64_t colorset_sampling_distance,
                        GGCAT_unitig_database& unitig_database){

        write_log("Building the color mapping", LogLevel::MAJOR);

        // Set up the worker threads that process the colored unitigs
        Sparse_Uint_Array_Builder builder(SBWT.number_of_subsets(), ram_bytes, n_threads);
        typename UnitigWorker::context_t context = {&SBWT, colorset_sampling_distance, &builder, {}};
        vector<unique_ptr<UnitigWorker>> workers;
        vector<UnitigWorker*> worker_ptrs;
        for(int64_t i = 0; i < n_threads; i++){
            workers.push_back(make_unique<UnitigWorker>(context));
            worker_ptrs.push_back(workers.back().get());
        }

        ThreadPool<UnitigWorker,UnitigWorkBatch> TP(worker_ptrs, n_threads * 2); 

        iterate_unitigs(unitig_database, coloring, TP);

        coloring.index_ptr = &SBWT;
        coloring.node_id_to_color_set_id = builder.finish();
        coloring.sets.prepare_for_queries();        
    }

};