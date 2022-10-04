#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <set>
#include <algorithm>
#include <filesystem>
#include "globals.hh"
#include "new_coloring.hh"
#include "WorkDispatcher.hh"
#include "EM_sort.hh"
#include "KMC_wrapper.hh"

using namespace std;

typedef int64_t LL; // long long

template <typename T>
set<T> intersect(const set<T>& S1, const set<T>& S2){
    set<T> ans;
    for(auto x : S1){
        if(S2.count(x) > 0){
            ans.insert(x);
        }
    }
    return ans;
}

// Stores the intersection into buf1 and returns the number of elements in the
// intersection (does not resize buf1). Buffer elements must be sorted.
// Assumes all elements in a buffer are distinct
LL intersect_buffers(vector<LL>& buf1, LL buf1_len, vector<LL>& buf2, LL buf2_len);

// Stores the union into buf1 and returns the number of elements in the
// intersection (does not resize buf1). Buffers elements must be sorted.
// Assumes all elements in a buffer are distinct. result_buf must have enough
// space to accommodate the union
LL union_buffers(vector<LL>& buf1, LL buf1_len, vector<LL>& buf2, LL buf2_len, vector<LL>& result_buf);

class Themisto{

public:

    class AlignerThread : public DispatcherConsumerCallback{

    private:

        AlignerThread(const AlignerThread&); // No copying
        AlignerThread& operator=(const AlignerThread&); // No copying

        // result must have at least the same length as S.
        LL get_rc(const char* S, LL S_size, string& result){
            for(LL i = 0; i < S_size; i++){
                result[S_size-1-i] = get_rc(S[i]);
            }
            return S_size;
        }

        char get_rc(char c){
            switch(c){
                case 'A': return 'T';
                case 'T': return 'A';
                case 'C': return 'G';
                case 'G': return 'C';
                default: return c;
            }
        }

    public:

        Themisto* kl;
        ParallelBaseWriter* out;
        bool reverse_complements;
        LL output_buffer_size;
        string output_buffer; // For printing
        vector<LL> temp_colorset_id_buffer;
        vector<LL> temp_colorset_id_buffer_rc;
        vector<LL> colorset_buffer;
        vector<LL> colorset_rc_buffer;
        vector<LL> temp_buffer;
        vector<LL> union_buffer;
        string rc_buffer;
        LL k;

        AlignerThread(Themisto* kl, ParallelBaseWriter* out, bool reverse_complements, LL output_buffer_size){
            this->kl = kl;
            this->out = out;
            this->reverse_complements = reverse_complements;
            this->output_buffer_size = output_buffer_size;
            this->colorset_buffer.resize(kl->coloring.n_colors);
            this->colorset_rc_buffer.resize(kl->coloring.n_colors);
            this->temp_buffer.resize(kl->coloring.n_colors);
            this->union_buffer.resize(kl->coloring.n_colors);
            this->k = kl->boss.get_k();
        }

        // Returns number of elements placed into putput
        LL intersect_colorsets(vector<LL>& ids, LL query_size, vector<LL>& output, vector<LL>& temp){
            bool found_at_least_one_nonempty_colorset = false;
            LL output_buf_size = 0;
            for(LL i = 0; i < query_size - k + 1; i++){
                if(ids[i] != -1 && (i == 0 || ids[i-1] != ids[i])){
                    // Nonempty and nonredundant color set
                    LL temp_buf_size = kl->coloring.get_colorset_to_buffer_by_id(ids[i], temp);
                    if(!found_at_least_one_nonempty_colorset){
                        // First color set that was found. Put it into output_buf
                        for(LL j = 0; j < temp_buf_size; j++) output[j] = temp[j];
                        output_buf_size = temp_buf_size;
                    } else{
                        // Intersect the colors in temp_buf with output_buf
                        output_buf_size = intersect_buffers(output, output_buf_size, temp, temp_buf_size);
                        if(output_buf_size == 0) return 0;
                    }
                    found_at_least_one_nonempty_colorset = true;
                }
            }
            return output_buf_size;
        }

        // Assumes temp_colorset_id_buffer and temp_colorset_id_buffer_rc are populated
        // Returns interection size
        LL do_intersections_with_legacy_behaviour(LL query_size){
            if(!reverse_complements){
                return intersect_colorsets(temp_colorset_id_buffer, query_size, colorset_buffer, temp_buffer);
            } else{
                LL forward_size = intersect_colorsets(temp_colorset_id_buffer, query_size, colorset_buffer, temp_buffer);
                LL rc_size = intersect_colorsets(temp_colorset_id_buffer_rc, query_size, colorset_rc_buffer, temp_buffer);
                LL union_size = union_buffers(colorset_buffer, forward_size,
                                              colorset_rc_buffer, rc_size, 
                                              union_buffer);
                for(LL i = 0; i < union_size; i++) colorset_buffer[i] = union_buffer[i];
                return union_size;
            }    
        }

        // Assumes temp_colorset_id_buffer is populated (and temp_colorset_id_buffer_rc if reverse complements
        // are enabled). Puts final intersected color set into colorset_buffer. Returns the number of elements
        // placed into colorset_buffer.
        LL do_intersections(LL query_size){
            if(!reverse_complements) {
                return intersect_colorsets(temp_colorset_id_buffer, query_size, colorset_buffer, temp_buffer);
            }

            // reverse complements are enbaled
            bool found_at_least_one_nonempty_colorset = false;
            LL temp_buf_size = 0;
            for(LL i = 0; i < query_size - k + 1; i++){
                bool need_to_intersect = true;
                if(temp_colorset_id_buffer[i] == -1 && temp_colorset_id_buffer_rc[query_size - k - i] == -1)
                    need_to_intersect = false; // Empty set
                else if(i > 0 && temp_colorset_id_buffer[i] == temp_colorset_id_buffer[i-1]
                         && temp_colorset_id_buffer_rc[query_size - k - i] == temp_colorset_id_buffer_rc[query_size - k - i + 1]){
                    need_to_intersect = false; // Same color set as previous
                }
                if(need_to_intersect){
                    LL forward_size = kl->coloring.get_colorset_to_buffer_by_id(temp_colorset_id_buffer[i], colorset_buffer);
                    LL rc_size = kl->coloring.get_colorset_to_buffer_by_id(temp_colorset_id_buffer_rc[query_size - k - i], colorset_rc_buffer);

                    LL union_size = union_buffers(colorset_buffer, forward_size,
                                                colorset_rc_buffer, rc_size, 
                                                union_buffer);
                    
                    if(!found_at_least_one_nonempty_colorset){
                        // First color set that was found. Put it into temp buf
                        for(LL j = 0; j < union_size; j++) temp_buffer[j] = union_buffer[j];
                        temp_buf_size = union_size;
                    } else{
                        // Intersect the colors in temp_buf with the union buffer
                        temp_buf_size = intersect_buffers(temp_buffer, temp_buf_size, union_buffer, union_size);
                        if(temp_buf_size == 0) return 0;
                    }
                    found_at_least_one_nonempty_colorset = true;
                }
            }

            // Copy to output
            for(LL i = 0; i < temp_buf_size; i++) colorset_buffer[i] = temp_buffer[i];
            
            return temp_buf_size;

        }


        virtual void callback(const char* S, LL S_size, int64_t string_id){
            if(S_size < k){
                write_log("Warning: query is shorter than k", LogLevel::MINOR);
                output_buffer += std::to_string(string_id) + "\n";
            }
            else{
                // Resize the buffers if they are too small
                if(this->temp_colorset_id_buffer.size() < S_size - k + 1){
                    temp_colorset_id_buffer.resize(S_size - k + 1);
                    if(reverse_complements) temp_colorset_id_buffer_rc.resize(S_size - k + 1);
                }

                // Get the color sets
                kl->get_nonempty_colorset_ids(S,S_size,temp_colorset_id_buffer);
                if(reverse_complements){
                    if(rc_buffer.size() < S_size) rc_buffer.resize(S_size);
                    get_rc(S, S_size, rc_buffer);
                    kl->get_nonempty_colorset_ids(rc_buffer.c_str(), S_size,temp_colorset_id_buffer_rc);
                }
                LL colorset_size = do_intersections(S_size);
                output_buffer += std::to_string(string_id) + " ";
                for(LL i = 0; i < colorset_size; i++)
                    output_buffer += std::to_string(colorset_buffer[i]) + " ";
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

    BOSS<sdsl::bit_vector> boss;
    Coloring coloring;
    vector<string> colors; // sequence id -> color name

    void construct_boss(string fastafile, LL k, LL memory_bytes, LL n_threads, bool revcomps){
/*
        write_log("Building KMC database", LogLevel::MAJOR);
        //string KMC_db_path_prefix = get_temp_file_manager().create_filename("KMC-");
        string KMC_db_path_prefix; LL n_kmers;
        std::tie(KMC_db_path_prefix, n_kmers) = run_kmc({fastafile}, k+1, n_threads, memory_bytes / (1LL << 30), 1, 1000000000);
        //KMC_wrapper(k+1, max(1LL, memory_bytes / (1LL << 30)), n_threads, fastafile, get_temp_file_manager().get_dir(), KMC_db_path_prefix, revcomps, get_log_level() == LogLevel::OFF);
        write_log("Building KMC database finished", LogLevel::MAJOR);
        Kmer_stream_from_KMC_DB edgemer_stream(KMC_db_path_prefix, revcomps);
        BOSS_builder<BOSS<sdsl::bit_vector>, Kmer_stream_from_KMC_DB> builder;
        write_log("Building BOSS from KMC database", LogLevel::MAJOR);
        boss = builder.build(edgemer_stream, memory_bytes, n_threads);

        // Delete the KMC database files. The temp file manager can not do this because
        // KMC appends suffixes to the filename and the manager does not know about that.
        write_log("Deleting KMC database", LogLevel::MAJOR);
        
        std::filesystem::remove(KMC_db_path_prefix + ".kmc_pre");
        std::filesystem::remove(KMC_db_path_prefix + ".kmc_suf");
*/
    }

    // Assumes boss has been built
    // colorfile: A filename for a file with one line for each sequence in fastafile. The line has the name of the color (a string).
    // If colorfile is emptry string, generates colors automatically
    void construct_colors(string fastafile, string colorfile, LL ram_bytes, LL n_threads, LL colorset_sampling_distance){
        assert(colorfile != "");
        vector<LL> seq_to_color = read_colorfile(colorfile);
        coloring.add_colors(boss, fastafile, seq_to_color, ram_bytes, n_threads, colorset_sampling_distance);
    }

    void save_boss(string filename){
        throwing_ofstream out(filename, ios::binary);
        boss.serialize(out.stream);
    }

    void save_colors(string filename){
        throwing_ofstream out(filename, ios::binary);
        coloring.serialize(out.stream);
    }

    void load_boss(string filename){
        throwing_ifstream in(filename, ios::binary);
        boss.load(in.stream);
    }

    void load_colors(string filename){
        throwing_ifstream in(filename, ios::binary);
        coloring.load(in.stream);
    }

    void save(string filename_prefix){
        save_boss(filename_prefix + ".tdbg");
        save_colors(filename_prefix + ".tcolors");
    }

    void load(string filename_prefix){
        load_boss(filename_prefix + ".tdbg");
        load_colors(filename_prefix + ".tcolors");
    }

    string to_string(){
        string boss_string = boss.to_string();
        string coloring_string = coloring.to_string(boss);
        return "--boss--\n" + boss_string + "\n--coloring--\n" + coloring_string;
    }

    // Returns start index, node_id
    // If not found, returns (-1, -1)
    pair<LL,LL> find_first_matching_kmer(const char* S, LL starting_from, LL end, LL k){
        for(LL i = starting_from; i <= end; i++){
            if(i + k - 1 <= end){
                LL node = boss.find_kmer(S + i);
                if(node != -1) return {i,node};
            }
        }
        return {-1,-1};
    }

    // Legacy function. Super slow because allocates new buffers
    // every time. Use pseudoalign_to_buffer instead
    set<LL> pseudoalign(const string& Q){
        vector<LL> output_buf(coloring.n_colors);
        vector<LL> temp_buf(coloring.n_colors);
        vector<LL> temp_buf2(Q.size());
        LL output_size = pseudoalign_to_buffer(Q.c_str(), Q.size(), temp_buf2, temp_buf, output_buf);
        set<LL> ans(output_buf.begin(), output_buf.begin() + output_size);
        return ans;
    }

    // Output buf should at least Q_size - k + 1 elements. Fills output buf so that element i
    // is the id of the color set of the i-th k-mer of Q. If the i-th k-mer does not exist in
    // the index, then the color set id is defined to be -1.
    void get_nonempty_colorset_ids(const char* Q, LL Q_size, vector<LL>& output_buf){

        LL k = boss.get_k();

        // Initialize outpuf buf values to -1
        for(LL i = 0; i < Q_size - k + 1; i++) output_buf[i] = -1;

        LL idx, node;
        std::tie(idx,node) = find_first_matching_kmer(Q, 0, Q_size-1, k);

        if(idx == -1) return; // None of the k-mers match the index
        
        output_buf[idx] = coloring.get_colorset_id(node, boss);

        while(idx + k - 1 < Q_size){
            // Loop invariant here: Q[idx..(idx+k-1)] is found in the index.

            if(idx + k >= Q_size) break; // At the last kmer of the query

            node = boss.walk(node, Q[idx+k]); // Try to walk to the next kmer
            if(node != -1){ // Success
                idx++;
                if(coloring.is_redundant(node) && output_buf[idx-1] != -1){
                    output_buf[idx] = output_buf[idx-1];
                } else {
                    output_buf[idx] = coloring.get_colorset_id(node, boss);
                }
                
            } else{ // Could not walk to the next kmer
                std::tie(idx,node) = find_first_matching_kmer(Q, idx+1, Q_size-1, k);
                if(idx == -1) break; // No more matching k-mers in this query
                output_buf[idx] = coloring.get_colorset_id(node, boss);
            }
        }
    }

    // Buffers must have enough space to accommodate all colors.
    // Puts the colors into the output buffer and returns the number
    // of elements. Does not resize the output buffer. colorset_id_buf should
    // have size at least Q_size - k + 1. output_buf and temp_buf should have at 
    // least size equal to the number of distinct colors.
    LL pseudoalign_to_buffer(const char* Q, LL Q_size, vector<LL>& colorset_id_buf, vector<LL>& temp_buf, vector<LL>& output_buf){

        LL k = boss.get_k();
        get_nonempty_colorset_ids(Q,Q_size,colorset_id_buf);

        bool found_at_least_one_nonempty_colorset = false;
        LL output_buf_size = 0;
        for(LL i = 0; i < Q_size - k + 1; i++){
            if(colorset_id_buf[i] != -1 && (i == 0 || colorset_id_buf[i-1] != colorset_id_buf[i])){
                // Nonempty and nonredundant color set
                LL temp_buf_size = coloring.get_colorset_to_buffer_by_id(colorset_id_buf[i], temp_buf);
                if(!found_at_least_one_nonempty_colorset){
                    // First color set that was found. Put it into output_buf
                    for(LL j = 0; j < temp_buf_size; j++) output_buf[j] = temp_buf[j];
                    output_buf_size = temp_buf_size;
                } else{
                    // Intersect the colors in temp_buf with output_buf
                    output_buf_size = intersect_buffers(output_buf, output_buf_size, temp_buf, temp_buf_size);
                    if(output_buf_size == 0) return 0;
                }
                found_at_least_one_nonempty_colorset = true;
            }
        }
        return output_buf_size;
    }    

    bool getline(throwing_ifstream& is, string& line){
        return is.getline(line);
    }

    bool getline(zstr::ifstream& is, string& line){
        if(std::getline(is,line)) return true;
        return false;
    }


    void pseudoalign_parallel(LL n_threads, Sequence_Reader_Buffered& sr, string outfile, bool reverse_complements, LL buffer_size, bool gzipped_output, bool sort_after){
        ParallelBaseWriter* out = nullptr;
	if (!outfile.empty()) {
	    if(gzipped_output) out = new ParallelGzipWriter(outfile);
	    else out = new ParallelOutputWriter(outfile);
	} else {
	    if(gzipped_output) out = new ParallelGzipWriter(cout);
	    else out = new ParallelOutputWriter(cout);
	}

        vector<DispatcherConsumerCallback*> threads;
        for(LL i = 0; i < n_threads; i++){
            AlignerThread* T = new AlignerThread(this, out, reverse_complements, buffer_size);
            threads.push_back(T);
        }

        run_dispatcher(threads, sr, buffer_size);

        // Clean up
        for(DispatcherConsumerCallback* t : threads) delete t;
        out->flush();
        delete out;

        if(sort_after){
            write_log("Sorting output file", LogLevel::MAJOR);
            string tempfile = get_temp_file_manager().create_filename("results_temp");
            if(gzipped_output){
                zstr::ifstream instream(outfile);
                zstr::ofstream outstream(tempfile);
                sort_parallel_output_file(instream, outstream);
            } else{
                throwing_ifstream instream(outfile);
                throwing_ofstream outstream(tempfile);
                sort_parallel_output_file(instream, outstream);
	    }
            std::filesystem::rename(tempfile, outfile);
        }
    }

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


    size_t get_file_size(const char* filename) {
        struct stat st;
        stat(filename, &st);
        return st.st_size;
    }

    // returns a vector where element i is the ref ids aligned with query i
    static vector<set<LL> > parse_output_format_from_disk(string filename){
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

};
