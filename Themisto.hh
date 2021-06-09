#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <set>
#include <algorithm>
#include "BD_BWT_index.hh"
#include "input_reading.hh"
#include "test_tools.hh"
#include "globals.hh"
#include "libwheeler/BOSS.hh"
#include "libwheeler/BOSS_builder.hh"
#include "Coloring.hh"
#include "NameMapping.hh"
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
LL intersect_buffers(vector<LL>& buf1, LL buf1_len, vector<LL>& buf2, LL buf2_len){

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

// Stores the union into buf1 and returns the number of elements in the
// intersection (does not resize buf1). Buffers elements must be sorted.
// Assumes all elements in a buffer are distinct. result_buf must have enough
// space to accommodate the union
LL union_buffers(vector<LL>& buf1, LL buf1_len, vector<LL>& buf2, LL buf2_len, vector<LL>& result_buf){

    auto end = std::set_union(
                    buf1.begin(), buf1.begin() + buf1_len, 
                    buf2.begin(), buf2.begin() + buf2_len, 
                    result_buf.begin()
    );
    return end - result_buf.begin();
}

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
            if(this->temp_colorset_id_buffer.size() < S_size - k + 1){
                temp_colorset_id_buffer.resize(S_size - k + 1);
                if(reverse_complements) temp_colorset_id_buffer_rc.resize(S_size - k + 1);
            }

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
    NameMapping name_mapping; // color name <-> color id
    vector<string> colors; // sequence id -> color name

    void construct_boss(string fastafile, LL k, LL memory_bytes, LL n_threads, bool revcomps){

        write_log("Building KMC database");
        string KMC_db_path_prefix = KMC_wrapper(k+1, max(1LL, memory_bytes / (1LL << 30)), n_threads, fastafile, temp_file_manager.get_dir(), revcomps);
        Kmer_stream_from_KMC_DB edgemer_stream(KMC_db_path_prefix, revcomps);
        BOSS_builder<BOSS<sdsl::bit_vector>, Kmer_stream_from_KMC_DB> builder;
        write_log("Building BOSS from KMC database");
        boss = builder.build(edgemer_stream, memory_bytes);

    }

    string generate_colorfile(string fastafile){
        string colorfile = temp_file_manager.get_temp_file_name("");
        throwing_ofstream out(colorfile);
        Sequence_Reader fr(fastafile, FASTA_MODE);
        LL seq_id = 0;
        while(!fr.done()){
            fr.get_next_query_stream().get_all();
            out << seq_id << "\n";
            seq_id++;
        }
        return colorfile;
    }

    // Assumes boss has been built
    // colorfile: A filename for a file with one line for each sequence in fastafile. The line has the name of the color (a string).
    // If colorfile is emptry string, generates colors automatically
    void construct_colors(string fastafile, string colorfile, LL ram_bytes, LL n_threads, LL colorset_sampling_distance){
        if(colorfile == "") colorfile = generate_colorfile(fastafile);
        colors = parse_color_name_vector(colorfile);
        name_mapping.build_from_namefile(colorfile);
        coloring.add_colors(boss, fastafile, name_mapping.names_to_ids(colors), ram_bytes, n_threads, colorset_sampling_distance);
    }

    void save_boss(string path_prefix){
        boss.save_to_disk(path_prefix);
    }

    void save_colors(string path_prefix){
        coloring.save_to_disk(path_prefix + "index-");
        name_mapping.save_to_disk(path_prefix + "mapping-");
        throwing_ofstream out(path_prefix + "names.txt");
        for(string name : colors) out << name << "\n"; 
    }

    void load_boss(string path_prefix){
        boss.load_from_disk(path_prefix);
    }

    void load_colors(string path_prefix){
        coloring.load_from_disk(path_prefix + "index-");
        name_mapping.load_from_disk(path_prefix + "mapping-");

        throwing_ifstream in(path_prefix + "names.txt");

        colors.clear();
        string name;
        while(in.getline(name)){
            colors.push_back(name);
        }
    }

    vector<string> parse_color_name_vector(string colorfile){
        throwing_ifstream colorstream(colorfile);
        vector<string> result;
        string color_name;
        while(colorstream.getline(color_name)){
            result.push_back(color_name);
        }
        return result;
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


    void pseudoalign_parallel(LL n_threads, Sequence_Reader& sr, string outfile, bool reverse_complements, LL buffer_size, bool gzipped_output, bool sort_after){
        ParallelBaseWriter* out = nullptr;
        if(gzipped_output) out = new ParallelGzipWriter(outfile);
        else out = new ParallelOutputWriter(outfile);

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
            write_log("Sorting output file");
            string tempfile = temp_file_manager.get_temp_file_name("results_temp");
            if(gzipped_output){
                zstr::ifstream instream(outfile);
                zstr::ofstream outstream(tempfile);
                sort_parallel_output_file(instream, outstream);
            } else{
                throwing_ifstream instream(outfile);
                throwing_ofstream outstream(tempfile);
                sort_parallel_output_file(instream, outstream);
            }
            copy_file(tempfile, outfile, 1024*1024);
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
    vector<set<LL> > parse_output_format_from_disk(string filename){
        vector<set<LL> > results;
        check_readable(filename);
        throwing_ifstream input(filename);
        string line;
        
        LL line_number = 0;
        while(input.getline(line)){
            vector<LL> tokens = parse_tokens<LL>(line);
            assert(tokens.size() >= 1);
            LL query_id = tokens[0];
            assert(query_id == line_number);
            set<LL> alignments;
            for(LL i = 1; i < tokens.size(); i++){
                alignments.insert(tokens[i]);
            }
            results.push_back(alignments);
            line_number++;
        }
        return results;
    }

};

class Themisto_Tester{

public:

    class TestCase{
        public:
        vector<string> genomes;
        
        unordered_map<string, set<string> > node_to_color_names; // kmer- > set of names
        vector<string> queries;
        vector<string> colex_kmers;
        LL n_colors;
        LL k;

        NameMapping nm;
        vector<string> colors; // for each sequence: the color name of that sequence
    };

    LL random_seed = 123674;


    string get_rc(string S){
        string R(S.rbegin(), S.rend());
        for(LL i = 0; i < R.size(); i++)
            R[i] = get_rc(R[i]);
        return R;
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

    vector<TestCase> generate_testcases(LL genome_length, LL n_genomes, LL n_queries, LL query_length, LL min_k, LL max_k, LL n_colors){
        srand(random_seed);
        vector<TestCase> testcases;

        for(LL rep = 0; rep < 5; rep++){
            for(LL k = min_k; k <= max_k; k++){
                TestCase tcase;
                tcase.k = k;

                vector<string> distinct_color_names;
                for(LL i = 0; i < n_colors; i++){
                    distinct_color_names.push_back(get_random_dna_string(10,2) + to_string(i));
                }

                // Build genomes and assign color ids
                for(LL i = 0; i < n_genomes; i++){
                    tcase.genomes.push_back(get_random_dna_string(genome_length,2));
                    tcase.colors.push_back(distinct_color_names[rand() % n_colors]);
                }

                tcase.nm.build_from_name_vector(tcase.colors);
                tcase.n_colors = tcase.nm.id_to_name.size();

                // Get all k-mers and colex-sort them
                set<string> all_kmers;
                for(LL i = 0; i < n_genomes; i++){
                    for(string kmer : get_all_distinct_kmers(tcase.genomes[i], k)) all_kmers.insert(kmer);
                }
                
                tcase.colex_kmers = vector<string>(all_kmers.begin(), all_kmers.end());
                sort(tcase.colex_kmers.begin(), tcase.colex_kmers.end(), colex_compare);

                // List k-mer sets for each color
                vector<set<string> > color_to_kmer_set(n_colors);
                for(LL genome_id = 0; genome_id < tcase.genomes.size(); genome_id++){
                    string genome = tcase.genomes[genome_id];
                    string colorname = tcase.colors[genome_id];
                    LL color_id = tcase.nm.name_to_id[colorname];
                    for(string kmer : get_all_distinct_kmers(genome,k)){
                        color_to_kmer_set[color_id].insert(kmer);
                    }
                }

                // List all color names for each k-mer (!= each node)
                for(string kmer : tcase.colex_kmers){
                    set<string> colorset;
                    for(LL color_id = 0; color_id < n_colors; color_id++){
                        if(color_to_kmer_set[color_id].count(kmer)){
                            colorset.insert(tcase.nm.id_to_name[color_id]);
                        }
                    }
                    tcase.node_to_color_names[kmer] = colorset;
                }
                

                // Build queries
                for(LL i = 0; i < n_queries; i++) tcase.queries.push_back(get_random_dna_string(query_length,2));

                testcases.push_back(tcase);
            }
        }
        return testcases;
    }

    void test_pseudoalign(string temp_dir){
        cerr << "Testing pseudolign" << endl;

        ColoringTester ct;
        LL testcase_id = -1;
        LL ref_length = 100;
        LL n_refs = 50;
        LL n_threads = 3;
        LL n_queries = 10000;
        LL query_length = 20;
        LL k_min = 1;
        LL k_max = 20;
        LL n_colors = 5;
        for(TestCase tcase : generate_testcases(ref_length,n_refs,n_queries,query_length,k_min,k_max,n_colors)){
            testcase_id++;
            cerr << "Running coloring testcase" << endl;
            
            ColoringTester::TestCase ct_testcase 
                = ct.generate_testcase(tcase.genomes, tcase.nm.names_to_ids(tcase.colors), tcase.k);
            ct.run_testcase(ct_testcase); // Check that this works

            cerr << "Running alignment testcase" << endl;
            throwing_ofstream genomes_out(temp_dir + "/genomes.fna");
            for(string genome : tcase.genomes){
                genomes_out << ">\n" << genome << "\n";
            }
            genomes_out.close();

            throwing_ofstream colors_out(temp_dir + "/colors.txt");
            for(LL i = 0; i < tcase.colors.size(); i++){
                colors_out << tcase.colors[i] << "\n";
            }
            colors_out.close();

            throwing_ofstream queries_out(temp_dir + "/queries.fna");
            for(string query : tcase.queries){
                queries_out << ">\n" << query << "\n";
            }
            queries_out.close();

            NameMapping nm(temp_dir + "/colors.txt");

            {
                Themisto kl_build;
                kl_build.construct_boss(temp_dir + "/genomes.fna", tcase.k, 1000, 2, true);
                kl_build.construct_colors(temp_dir + "/genomes.fna", temp_dir + "/colors.txt", 1000, 3, 10);
                kl_build.save_boss(temp_dir + "/boss-");
                kl_build.save_colors(temp_dir + "/colors-");
            }

            Themisto kl; // Load back to test serialization
            kl.load_boss(temp_dir + "/boss-");
            kl.load_colors(temp_dir + "/colors-");

            // Check against reference boss
            BOSS<sdsl::bit_vector> ref_boss = build_BOSS_with_maps(tcase.genomes, tcase.k, true); // Reverse complements enabled because the KMC construction uses those
            assert(ref_boss == kl.boss);

            // Run without rc
            string final_file = temp_file_manager.get_temp_file_name("finalfile");
            Sequence_Reader sr(temp_dir + "/queries.fna", FASTA_MODE);
            kl.pseudoalign_parallel(n_threads, sr, final_file, false, 300, false, true);
            vector<set<LL> > our_results = kl.parse_output_format_from_disk(final_file);

            // Run with rc
            string final_file_rc = temp_file_manager.get_temp_file_name("finalfile_rc");
            Sequence_Reader sr2(temp_dir + "/queries.fna", FASTA_MODE);
            kl.pseudoalign_parallel(n_threads, sr2, final_file_rc, true, 300, false, true);
            vector<set<LL> > our_results_rc = kl.parse_output_format_from_disk(final_file_rc);

            //cout << our_results << endl;

            for(LL i = 0; i < tcase.queries.size(); i++){
                string query = tcase.queries[i];
                set<string> brute = pseudoalign_to_colors_trivial(query, tcase, false);
                set<string> brute_rc = pseudoalign_to_colors_trivial(query, tcase, true);
                
                set<string> nonbrute;
                for(LL id : our_results[i]) nonbrute.insert(nm.id_to_name[id]);

                set<string> nonbrute_rc;
                for(LL id : our_results_rc[i]) nonbrute_rc.insert(nm.id_to_name[id]);

                assert(brute == nonbrute);
                assert(brute_rc == nonbrute_rc);
            }
        }

    }

    // Returns set of color names
    set<string> pseudoalign_to_colors_trivial(string& query, TestCase& tcase, bool reverse_complements){
        set<string> alignments;
        for(LL i = 0; i < tcase.genomes.size(); i++) alignments.insert(tcase.colors[i]); // All color names

        bool at_least_one = false;
        // For each k-mer in query, get the color name set and intersect that with alignments
        for(string kmer : get_all_kmers(query, tcase.k)){
            set<string> names = tcase.node_to_color_names[kmer];
            if(reverse_complements){
                set<string> names_rc = tcase.node_to_color_names[get_rc(kmer)];
                for(string name : names_rc) names.insert(name);
            }
            if(names.size() >= 1) {
                at_least_one = true;
                alignments = intersect(alignments, names);
            }
        }

        if(at_least_one == false) alignments.clear();
        return alignments;
    }

};

void test_pseudoalign(string temp_dir){
    Themisto_Tester tester;
    tester.test_pseudoalign(temp_dir);
}
