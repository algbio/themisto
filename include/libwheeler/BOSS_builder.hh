#pragma once

#include "Kmer.hh"
#include "BOSS.hh"
#include "EM_sort.hh"
#include "stdafx.h"
#include "kmc_file.h"
#include "test_tools.hh"
#include <unordered_map>
#include <stdexcept>

template<typename iter_t>
class Kmer_stream_from_iterator_pair{

public:

    iter_t it;
    iter_t end;

    Kmer_stream_from_iterator_pair(iter_t begin, iter_t end) : it(begin), end(end){}

    bool done(){
        return it == end;
    }

    Kmer<KMER_MAX_LENGTH> next(){
        Kmer<KMER_MAX_LENGTH> x = *it;
        it++;
        return x;
    }

};

class Kmer_stream_in_memory{
public:

    const vector<string>& reads;
    LL cur_read_id = 0;
    LL cur_read_pos = 0;
    LL k;

    /**
     * Assumers reads have characters only from alphabet {A,C,G,T}
     */
    Kmer_stream_in_memory(const vector<string>& reads, LL k) : reads(reads), k(k){
        
        cur_read_id = reads.size();

        // Find the first read that has length at least k
        for(LL i = 0; i < (LL)reads.size(); i++){
            if((LL)reads[i].size() >= k){
                cur_read_id = i;
                break;
            }
        }
    }

    bool done(){
        return cur_read_id == (LL)reads.size();
    }

    Kmer<KMER_MAX_LENGTH> next(){
        Kmer<KMER_MAX_LENGTH> x(reads[cur_read_id].substr(cur_read_pos++, k));

        // Get ready for the next read
        while(cur_read_id < (LL)reads.size() && cur_read_pos + k - 1 >= (LL)reads[cur_read_id].size()){
            cur_read_id++;
            cur_read_pos = 0;
        }

        return x;
    }
};

// Also gives reverse compelments
class Kmer_stream_from_KMC_DB{

private:

CKMCFile kmer_database;
CKmerAPI kmer_object;

uint32 _kmer_length;
uint32 _mode;
uint32 _counter_size;
uint32 _lut_prefix_length;
uint32 _signature_len;
uint32 _min_count;
uint64 _max_count;
uint64 _total_kmers;

bool add_revcomps;
std::string str;
std::string str_revcomp;
bool revcomp_next = false;

char get_rc(char c){
    switch(c){
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: cerr << "Error getting reverse complement from " << c << endl; exit(1);
    }   
}

void reverse_complement(string& S){
    std::reverse(S.begin(), S.end());
    for(char& c : S) c = get_rc(c);
}   

public:

    Kmer_stream_from_KMC_DB(string KMC_db_path, bool add_revcomps) : add_revcomps(add_revcomps) {
        if (!kmer_database.OpenForListing(KMC_db_path)){
            throw std::runtime_error("Error opening KMC database " + KMC_db_path);
        }

		kmer_database.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count, _max_count, _total_kmers);

        kmer_object = CKmerAPI(_kmer_length);
	
    }

    bool done(){
        return (!add_revcomps || !revcomp_next) && kmer_database.Eof();
    }

    Kmer<KMER_MAX_LENGTH> next(){
        if(add_revcomps && revcomp_next){
            revcomp_next = false;
            return Kmer<KMER_MAX_LENGTH>(str_revcomp);
        }

        float counter_f;
        uint32 counter_i;
		if(_mode){ //quake compatible mode
			kmer_database.ReadNextKmer(kmer_object, counter_f);
		}
		else {			
			kmer_database.ReadNextKmer(kmer_object, counter_i);
		}

        kmer_object.to_string(str);
        if(add_revcomps){
            str_revcomp = str;
            reverse_complement(str_revcomp);
            if(str != str_revcomp) revcomp_next = true;
        }

        return Kmer<KMER_MAX_LENGTH>(str);

    }

};

template<typename payload_t>
class Kmer_sorter_disk{

private:

// Have pointer members -> no copying
Kmer_sorter_disk(Kmer_sorter_disk const&) = delete;
Kmer_sorter_disk& operator=(Kmer_sorter_disk const&) = delete;

bool is_sorted = false;
bool all_read = false;
string unsorted_filename;
string sorted_filename;
Buffered_ofstream unsorted;
Buffered_ifstream sorted;
LL mem_budget_bytes = 1e9;
LL n_threads;
char* out_buffer;
char* in_buffer;

pair<Kmer<KMER_MAX_LENGTH>, payload_t> top;

void update_top(){
    Kmer<KMER_MAX_LENGTH> x; payload_t E;

    sorted.read(in_buffer, Kmer<KMER_MAX_LENGTH>::size_in_bytes());
    x.load(in_buffer);

    if(sorted.eof()){
        all_read = true;
        return;
    }

    sorted.read(in_buffer, payload_t::size_in_bytes());
    E.load(in_buffer);

    top = {x,E};
}

public:

    Kmer_sorter_disk(LL n_threads) : n_threads(n_threads) {
        unsorted_filename = get_temp_file_manager().create_filename("kmers");
        sorted_filename = get_temp_file_manager().create_filename("kmers_sorted");
        unsorted.open(unsorted_filename, ios_base::binary);
        LL buffer_size = max(Kmer<KMER_MAX_LENGTH>::size_in_bytes(), payload_t::size_in_bytes());
        out_buffer = (char*)malloc(buffer_size);
        in_buffer = (char*)malloc(buffer_size);
    }

    void set_mem_budget(LL bytes){
        mem_budget_bytes = bytes;
    }

    void add(Kmer<KMER_MAX_LENGTH> x, payload_t E){

        x.serialize(out_buffer);
        unsorted.write(out_buffer, Kmer<KMER_MAX_LENGTH>::size_in_bytes());

        E.serialize(out_buffer); // Overwrites x
        unsorted.write(out_buffer, payload_t::size_in_bytes());
    }

    void sort(){
        unsorted.close();
        auto cmp = [](const char* x, const char* y){
            Kmer<KMER_MAX_LENGTH> A, B;
            A.load(x); // x is the concatenation of the serialization of a Kmer and an edgeset
            B.load(y); // y is the concatenation of the serialization of a Kmer and an edgeset
            return A < B; // Colexicographic compare
        };
        EM_sort_constant_binary(unsorted_filename, sorted_filename, cmp, mem_budget_bytes, Kmer<KMER_MAX_LENGTH>::size_in_bytes() + payload_t::size_in_bytes(), n_threads);
        is_sorted = true;
        sorted.open(sorted_filename, ios_base::binary);
        update_top();
    }

    void reset_stream(){
        all_read = false;
        sorted.close();
        sorted.open(sorted_filename, ios_base::binary);
        update_top();
    }

    bool stream_done() const{
        return all_read;
    }

    pair<Kmer<KMER_MAX_LENGTH>, payload_t> stream_next(){
        assert(is_sorted);
        pair<Kmer<KMER_MAX_LENGTH>, payload_t> ret = top;
        update_top();
        return ret;
    }

    pair<Kmer<KMER_MAX_LENGTH>, payload_t> peek_next(){
        return top;
    }

    ~Kmer_sorter_disk(){
        free(out_buffer);
        free(in_buffer);
    }

};

template<typename payload_t>
class Kmer_sorter_in_memory{

private:

vector<pair<Kmer<KMER_MAX_LENGTH>, payload_t> > vec;
bool sorted = false;
LL stream_idx = 0;

public:

    void add(Kmer<KMER_MAX_LENGTH> x, payload_t E){
        vec.push_back({x,E});
    }

    void sort(){
        std::sort(vec.begin(), vec.end(), [](const pair<Kmer<KMER_MAX_LENGTH>, payload_t>& A, const pair<Kmer<KMER_MAX_LENGTH>, payload_t>& B){
            return A.first < B.first; // Colexicographic comparison
        });
        sorted = true;
    }

    // This does nothing but it's here to have the same interface as the other sorters.
    void set_mem_budget(LL bytes){
        (void)bytes;
        return; 
    }

    void reset_stream(){
        stream_idx = 0;
    }

    bool stream_done() const{
        return stream_idx == (LL)vec.size();
    }

    pair<Kmer<KMER_MAX_LENGTH>, payload_t> stream_next(){
        assert(sorted);
        return vec[stream_idx++];
    }

    pair<Kmer<KMER_MAX_LENGTH>, payload_t> peek_next(){
        assert(sorted);
        return vec[stream_idx];
    }


};

// A single sorted stream out of two sorted streams
template<typename payload_t, typename sorter_t>
class Kmer_stream_merger{

    sorter_t& A;
    sorter_t& B;

public:

    Kmer_stream_merger(sorter_t& A, sorter_t& B) : A(A), B(B){}

    bool stream_done(){
        return A.stream_done() && B.stream_done();
    }

    pair<Kmer<KMER_MAX_LENGTH>, payload_t> stream_next(){
        if(A.stream_done()) return B.stream_next();
        if(B.stream_done()) return A.stream_next();
        if(A.peek_next().first < B.peek_next().first) return A.stream_next();
        else return B.stream_next();
    }

};

// Edgemer stream should have functions:
// Kmer next(); // Returns the next (k+1)-mer
// bool done(); // Returns whether all (k+1)-mers have been iterated
template <typename boss_t, typename edgemer_stream>
class BOSS_builder{
    
private:

    struct Edgeset{
        //std::bitset<8> edges;
        uint8_t data;

        Edgeset() : data(0) {}

        void set_bit(uint8_t bit, uint8_t value){
            data = (data & ~(((uint8_t)1) << bit)) | (value << bit);
        }

        bool get_bit(uint8_t bit) const{
            return (data >> bit) & ((uint8_t)1);
        }

        bool have_out(char c) const{
            switch(c){
                case 'A': return get_bit(0);
                case 'C': return get_bit(1);
                case 'G': return get_bit(2);
                case 'T': return get_bit(3);
            }
            assert(false); // Invalid character
            return false;
        }

        bool have_in(char c) const{
            switch(c){
                case 'A': return get_bit(4);
                case 'C': return get_bit(5);
                case 'G': return get_bit(6);
                case 'T': return get_bit(7);
            }
            assert(false); // Invalid character
            return false;
        }

        void set_have_out(char c, bool have){
            switch(c){
                case 'A': set_bit(0, have); return;
                case 'C': set_bit(1, have); return;
                case 'G': set_bit(2, have); return;
                case 'T': set_bit(3, have); return;
            }
            assert(false); // Invalid character
        }

        void set_have_in(char c, bool have){
            switch(c){
                case 'A': set_bit(4, have); return;
                case 'C': set_bit(5, have); return;
                case 'G': set_bit(6, have); return;
                case 'T': set_bit(7, have); return;
            }
            assert(false); // Invalid character
        }

        string to_string(){
            stringstream ss;
            ss << "In: {";
            for(char c : string("ACGT")){
                if(have_in(c)) ss << c;
            }
            ss << "} Out: {";
            for(char c : string("ACGT")){
                if(have_out(c)) ss << c;
            }
            ss << "}";
            return ss.str();
        }

        static LL size_in_bytes(){
            return sizeof(data);
        }

        void serialize(ostream& out) const{
            out.write((char*)(&data), sizeof(data));
        }

        void load(istream& in){
            in.read((char*)(&data), sizeof(data));
        }

        // out must have at least size_in_bytes() bytes of space
        void serialize(char* out) const{
            memcpy(out, (char*)(&data), sizeof(data));
        }

        void load(const char* in){
            memcpy((char*)(&data), in, sizeof(data));
        }

    };

    // Returns number of prefixes added
    LL process_prefixes(Kmer<KMER_MAX_LENGTH> x, Edgeset E, Kmer_sorter_disk<Edgeset>& new_sorter){
        LL k = x.get_k();
        LL count = 0;
        if(!E.have_in('A') && !E.have_in('C') && !E.have_in('G') && !E.have_in('T')){
            // Add prefixes
            char c = '\0';
            while(true){
                Edgeset E_new;
                if(x.get_k() > 0)
                    E_new.set_have_in('A', true); // This is actually "dollar" but we use 'A' as dollar because the edgeset does not have a bit for dollar
                if(x.get_k() < k)
                    E_new.set_have_out(c, true);
                new_sorter.add(x, E_new);
                count++;
                if(x.get_k() > 0){
                    c = x.last();
                    x.dropright();
                } else break;
            }
        }
        return count;
    }

    // Assumes the old sorter is already sorted and the other is just initialized.
    // Adds the new edges to the new sorter. Returns number of dummies added.
    LL add_dummies(Kmer_sorter_disk<Edgeset>& old_sorter, Kmer_sorter_disk<Edgeset>& new_sorter){
        string ACGT = "ACGT";
        Kmer<KMER_MAX_LENGTH> prev_kmer;
        Edgeset cur_edgeset;
        LL loopcount = 0;
        LL addcount = 0;

        while(!old_sorter.stream_done()){
            Kmer<KMER_MAX_LENGTH> cur_kmer; Edgeset E;
            tie(cur_kmer,E) = old_sorter.stream_next();

            if(loopcount > 0 && cur_kmer != prev_kmer){
                // New node -> process the collected edges of the previous node
                addcount += process_prefixes(prev_kmer, cur_edgeset, new_sorter);
                cur_edgeset = Edgeset();
            }

            // Add the edge to collection
            for(char c : ACGT){
                if(E.have_in(c)) cur_edgeset.set_have_in(c, true);
                if(E.have_out(c)) cur_edgeset.set_have_out(c, true);
            }

            // Special case for processing the last node
            if(old_sorter.stream_done()) 
                addcount += process_prefixes(cur_kmer, cur_edgeset, new_sorter);

            prev_kmer = cur_kmer;
            loopcount++;
        }
        return addcount;
    }

public:

    boss_t build(edgemer_stream& input, LL mem_bytes, LL n_threads){
        if(input.done()) return boss_t(); 

        Kmer_sorter_disk<Edgeset> sorter1(n_threads);
        sorter1.set_mem_budget(mem_bytes);
        LL k = 0;
        LL n_records_written = 0;
        LL edge_count = 0;
        while(!input.done()){
            Kmer<KMER_MAX_LENGTH> edgemer = input.next();
            k = edgemer.get_k() - 1; // Edgemers are (k+1)-mers

            char first = edgemer.get(0);
            char last = edgemer.get(edgemer.get_k() - 1);

            Kmer<KMER_MAX_LENGTH> prefix = edgemer.copy().dropright();
            Edgeset prefixE; prefixE.set_have_out(last, true);
            sorter1.add(prefix, prefixE); n_records_written++;

            Kmer<KMER_MAX_LENGTH> suffix = edgemer.copy().dropleft();
            Edgeset suffixE; suffixE.set_have_in(first, true);
            sorter1.add(suffix, suffixE); n_records_written++;

            edge_count++;

        }
        write_log("Sorting " + to_string(edge_count) + " (k+1)-mers", LogLevel::MAJOR);
        sorter1.sort();

        // Dummies
        write_log("Adding dummy (k+1)-mers", LogLevel::MAJOR);
        Kmer_sorter_disk<Edgeset> sorter2(n_threads);
        sorter2.set_mem_budget(mem_bytes);
        LL n_dummies_disk = add_dummies(sorter1, sorter2);
        sorter2.sort();
        sorter1.reset_stream(); // Rewind back to start because add_dummies iterates over this

        // Build BOSS vectors
        Kmer_stream_merger<Edgeset, Kmer_sorter_disk<Edgeset>> merger(sorter1, sorter2);
        
        vector<bool> I, O;
        string outlabels;
        string ACGT = "ACGT";
        Kmer<KMER_MAX_LENGTH> prev_kmer;
        Edgeset cur_edgeset;
        LL loopcount = 0;

        if(n_dummies_disk == 0){
            // There were no dummies. But the BOSS structure wants there always
            // the be a source node. So we need to add it separately.
            I.push_back(1); // No in-edges
            O.push_back(1); // No out-edges
        }

        auto add_node = [&](){
            I.push_back(1);
            O.push_back(1);
            for(char c : ACGT){
                if(cur_edgeset.have_in(c)) I.push_back(0);
                if(cur_edgeset.have_out(c)){
                    O.push_back(0);
                    outlabels += c;
                }
            }
        };

        write_log("Constructing Wheeler BOSS components.", LogLevel::MAJOR);
        while(!merger.stream_done()){
            Kmer<KMER_MAX_LENGTH> cur_kmer; Edgeset E;
            tie(cur_kmer,E) = merger.stream_next();

            if(loopcount > 0 && cur_kmer != prev_kmer){
                // New node -> process the collected edges of a the previous node
                add_node();
                cur_edgeset = Edgeset();
            }

            // Add the edge to collection
            for(char c : ACGT){
                if(E.have_in(c)) cur_edgeset.set_have_in(c, true);
                if(E.have_out(c)){
                    cur_edgeset.set_have_out(c, true);
                }
            }

            // Special case for processing the last node
            if(merger.stream_done()) add_node();

            prev_kmer = cur_kmer;
            loopcount++;
        }

        assert(I.size() == O.size());

        return boss_t(outlabels, I, O, k);

    }

};

// Constructs the edge-centric de-bruijn graph where nodes are k-mers.
// i.e. for every (k+1)-mer, adds and edge between the prefix and suffix k-mer.
// The set of nodes in implicitly defined as the set of endpoints of all these
// nodes.
// Warning: this is very inefficiently implemented! Used for debug purposes
template<typename bitvector_t = sdsl::bit_vector>
BOSS<bitvector_t> build_BOSS_with_maps(vector<string> reads, LL k, bool include_reverse_complements){

    if(include_reverse_complements){
        LL n = reads.size();
        for(LL i = 0; i < n; i++){
            reads.push_back(get_rc(reads[i]));
        }
    }

    struct Edge_Info{
        set<char> inlabels;
        set<char> outlabels;
    };

    struct colex_compare_object{
        bool operator()(const string& A, const string& B) const{
            return colex_compare(A,B);
        }
    };

    map<string, Edge_Info, colex_compare_object> M; // edgemer -> edge info

    // Add all edges
    for(string seq : reads){
        if(seq.size() >= k+1){
            for(string x : get_all_distinct_kmers(seq, k+1)){
                M[x.substr(0,k)].outlabels.insert(x.back());
                M[x.substr(1,k)].inlabels.insert(x[0]);
            }
        }
    }


    // Add dummy nodes
    map<string, Edge_Info, colex_compare_object> M_copy = M; // Take a copy so we don't edit M while we iterate it
    for(auto P : M_copy){
        string kmer = P.first;
        if(P.second.inlabels.size() == 0){
            // Need to add the prefixes
            for(LL len = 0; len <= k; len++){
                string prefix = kmer.substr(0,len);
                if(len != 0) M[prefix].inlabels.insert('$');
                if(len != k) M[prefix].outlabels.insert(kmer[len]);
            }
        }
    }

    // Make sure the root node exists
    if(M.find("") == M.end()) M[""] = Edge_Info();

    write_log("map build added " + to_string(M.size() - M_copy.size()) + " dummies (k = " + to_string(k) + ")", LogLevel::MAJOR);


    // Build the boss data structures    
    string outlabels; // Concatenation in outedge label sets in colexicographic order
    vector<bool> outdegrees; // Concatenation of unary representations. Indegree d -> 1 0^d
    vector<bool> indegrees; // Concatenation of unary representations. Outdegree d -> 1 0^d

    for(auto P : M){
        outdegrees.push_back(1);
        string debug_outlabels;
        for(char c : P.second.outlabels){
            outlabels.push_back(c);
            outdegrees.push_back(0);
        }

        indegrees.push_back(1);
        for(char c : P.second.inlabels){
            (void)c; // Unsed variable. Make compiler happy
            if(P.first != "") indegrees.push_back(0);
        }         
    }

    return BOSS<bitvector_t>(outlabels, indegrees, outdegrees, k);
}