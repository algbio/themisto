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
#include "mark_kmers.hh"
#include "kmer_tools.hh"
#include "EM_sort.hh"

using namespace std;

typedef int64_t LL; // long long

// list_kmers is defined in KMC_wrapper.cpp
extern void list_kmers(int64_t k, int64_t ram_gigas, int64_t n_threads, 
                       string fastafile, string outfile, string tempdir);

class Rank_String{

// Represents a string and supports rank queries

private:

sdsl::wt_huff<sdsl::bit_vector> wt;

public:

    Rank_String(const string& S) {
        sdsl::construct_im(wt, S.c_str(),1); // 1: file format is a sequence, not a serialized sdsl object
    }

    // Rank up to i, excluding i. The indexing of i starts from 0.
    LL rank(LL i, char c){
        return wt.rank(i, c);
    }

    // Return the position of the i-th c. The indexing of i starts from 1.
    LL select(LL i, char c){
        return wt.select(i,c);
    }

    LL size(){
        return wt.size();
    }

    char char_at(LL i){
        return wt[i];
    }

    void save_to_disk(string path_prefix){
        sdsl::store_to_file(wt, path_prefix + "wt");
    }

    void load_from_disk(string path_prefix){
        check_readable(path_prefix + "wt");
        ifstream input(path_prefix + "wt");
        wt.load(input);
    }

};

class BOSS{

// Based on the Wheeler graph representation of the De Bruijn graph, see:
// Gagie, T., Manzini, G., & Sirén, J. (2017). Wheeler graphs: A framework for BWT-based data structures. Theoretical computer science, 698, 67-78.
// Represents all k-mers of the CYCLIC version of the input string.
// Since it's cyclic, every node is reachable from every node, and every node has at least one successor and predecessor

private:

    Rank_String* outlabels; // Should never be null
    sdsl::bit_vector indegs;
    sdsl::bit_vector outdegs;
    vector<LL> C; // Array of length 256. C[c] is the number of edges with label strictly less than c in the graph
    
    sdsl::rank_support_v<1> indegs_rs; // rank support
    sdsl::select_support_mcl<1> indegs_ss1; // select support
    sdsl::select_support_mcl<0> indegs_ss0; // select support
    
    sdsl::rank_support_v<1> outdegs_rs; // rank support
    sdsl::select_support_mcl<1> outdegs_ss1; // select support
    sdsl::select_support_mcl<0> outdegs_ss0; // select support

    vector<char> alphabet; // Not serialized to disk but build from outlabels

    LL n_nodes;
    LL k;
    

    void set_supports(){
        this->indegs_rs.set_vector(&(this->indegs)); 
        this->indegs_ss1.set_vector(&(this->indegs));
        this->indegs_ss0.set_vector(&(this->indegs));
        this->outdegs_rs.set_vector(&(this->outdegs)); 
        this->outdegs_ss1.set_vector(&(this->outdegs));
        this->outdegs_ss0.set_vector(&(this->outdegs));
    }

public:

    void save_to_disk(string path_prefix){
        outlabels->save_to_disk(path_prefix + "outlabels-");
        sdsl::store_to_file(indegs, path_prefix + "indegs");
        sdsl::store_to_file(outdegs, path_prefix + "outdegs");

        sdsl::store_to_file(indegs_rs, path_prefix + "indegs_rs");
        sdsl::store_to_file(indegs_ss1, path_prefix + "indegs_ss1");
        sdsl::store_to_file(indegs_ss0, path_prefix + "indegs_ss0");

        sdsl::store_to_file(outdegs_rs, path_prefix + "outdegs_rs");
        sdsl::store_to_file(outdegs_ss1, path_prefix + "outdegs_ss1");
        sdsl::store_to_file(outdegs_ss0, path_prefix + "outdegs_ss0");

        check_writable(path_prefix + "C");
        check_writable(path_prefix + "constants");
        check_writable(path_prefix + "alphabet");

        ofstream C_out(path_prefix + "C");
        for(auto x : C) C_out << x << " ";
        C_out << endl;
        C_out.close();

        ofstream constants_out(path_prefix + "constants");
        constants_out << n_nodes << " " << k << endl;
        constants_out.close();
    }

    void load_from_disk(string path_prefix){
        outlabels->load_from_disk(path_prefix + "outlabels-");

        check_readable(path_prefix + "indegs");
        check_readable(path_prefix + "outdegs");
        check_readable(path_prefix + "indegs_rs");
        check_readable(path_prefix + "indegs_ss1");
        check_readable(path_prefix + "indegs_ss0");
        check_readable(path_prefix + "outdegs_rs");
        check_readable(path_prefix + "outdegs_ss1");
        check_readable(path_prefix + "outdegs_ss0");
        check_readable(path_prefix + "C");
        check_readable(path_prefix + "constants");

        {
            ifstream input(path_prefix + "indegs");
            indegs.load(input);
        }

        {
            ifstream input(path_prefix + "outdegs");
            outdegs.load(input);
        }

        {
            ifstream input(path_prefix + "indegs_rs");
            indegs_rs.load(input);
        }

        {
            ifstream input(path_prefix + "indegs_ss1");
            indegs_ss1.load(input);
        }

        {
            ifstream input(path_prefix + "indegs_ss0");
            indegs_ss0.load(input);
        }

        {
            ifstream input(path_prefix + "outdegs_rs");
            outdegs_rs.load(input);
        }

        {
            ifstream input(path_prefix + "outdegs_ss1");
            outdegs_ss1.load(input);
        }

        {
            ifstream input(path_prefix + "outdegs_ss0");
            outdegs_ss0.load(input);
        }

        {
            C.resize(256);
            ifstream input(path_prefix + "C");
            for(LL i = 0; i < 256; i++){
                if(!input.good()){
                    cerr << "Error eading file: " << path_prefix + "C" << endl;
                    exit(1);
                }
                LL x; input >> x; 
                C[i] = x;
            }
        }

        {
            ifstream input(path_prefix + "constants");
            input >> n_nodes >> k; // todo: error handling
        }

        alphabet.clear();
        set<char> alphabet_set;
        for(LL i = 0; i < outlabels->size(); i++) alphabet_set.insert(outlabels->char_at(i)); // Todo: heavy
        for(char c : alphabet_set) alphabet.push_back(c);

    }

    BOSS(const BOSS& other){
        if(&other != this){

            delete this->outlabels;
            this->outlabels = new Rank_String(*other.outlabels);

            this->indegs = other.indegs;
            this->outdegs = other.outdegs;
            this->C = other.C;

            this->indegs_rs = other.indegs_rs;
            this->indegs_ss1 = other.indegs_ss1;
            this->indegs_ss0 = other.indegs_ss0;

            this->outdegs_rs = other.outdegs_rs;
            this->outdegs_ss1 = other.outdegs_ss1;
            this->outdegs_ss0 = other.outdegs_ss0;

            this->alphabet = other.alphabet;
            this->n_nodes = other.n_nodes;
            this->k = other.k;
        }
    }

    BOSS& operator=(const BOSS& other){
        // This function is almost identical to the copy constructor so one would think I could just call the 
        // copy constructor from here. But it's not that straightforward. Stackoverflow says that I can do 
        // that but then I need a swap-method (and the standard library swap won't do for some reason). But 
        // the stackoverflow answer does not give an example of how the swap should be implemented! 
        // So I'm just going with the duplicate code.
        // https://stackoverflow.com/questions/1734628/copy-constructor-and-operator-overload-in-c-is-a-common-function-possible
        if(&other != this){

            delete this->outlabels;
            this->outlabels = new Rank_String(*other.outlabels);

            this->indegs = other.indegs;
            this->outdegs = other.outdegs;
            this->C = other.C;

            this->indegs_rs = other.indegs_rs;
            this->indegs_ss1 = other.indegs_ss1;
            this->indegs_ss0 = other.indegs_ss0;

            this->outdegs_rs = other.outdegs_rs;
            this->outdegs_ss1 = other.outdegs_ss1;
            this->outdegs_ss0 = other.outdegs_ss0;
            
            this->alphabet = other.alphabet;
            this->n_nodes = other.n_nodes;
            this->k = other.k;
        }

        return *this;
    }

    BOSS(){
        string S = "";
        outlabels = new Rank_String(S);
    }

    BOSS(const string& outlabels_string, const sdsl::bit_vector& indegs, const sdsl::bit_vector& outdegs, const vector<LL>& C, LL k) 
    : indegs(indegs), outdegs(outdegs), C(C), k(k){
        outlabels = new Rank_String(outlabels_string);
        sdsl::util::init_support(indegs_rs, &indegs);
        sdsl::util::init_support(indegs_ss1, &indegs);
        sdsl::util::init_support(indegs_ss0, &indegs);
        sdsl::util::init_support(outdegs_rs, &outdegs);
        sdsl::util::init_support(outdegs_ss1, &outdegs);
        sdsl::util::init_support(outdegs_ss0, &outdegs);
        n_nodes = indegs_rs(indegs.size());
        set<char> alphabet_set(outlabels_string.begin(), outlabels_string.end());
        for(char c : alphabet_set) alphabet.push_back(c);
    }

    ~BOSS(){
        delete outlabels;
    }
    
    string get_outlabels(){
        string S;
        for(LL i = 0; i < outlabels->size(); i++) S += outlabels->char_at(i);
        return S;   
    }

    sdsl::bit_vector get_indegs(){return indegs;}
    sdsl::bit_vector get_outdegs(){return outdegs;}
    vector<LL> get_C_array(){return C;}
    vector<char> get_alphabet(){return alphabet;}
    LL get_number_of_nodes(){return n_nodes;}
    LL get_number_of_edges(){return outlabels->size();}
    LL get_k(){return k;}

    LL search(const char* kmer){
        // Make sure the support structures point to the right vector (if the object is copied the addressed can change)
        set_supports(); 

        // [l,r] = current colexicographic range (zero-indexed)
        LL l = 0;
        LL r = n_nodes-1;

        for(LL i = 0; i < k; i++){
            LL c = kmer[i];

            LL start = outdegs_ss1(l+1) - l;
            LL end;
            if(r == n_nodes-1) end = outlabels->size()-1;
            else end = outdegs_ss1(r+1+1) - (r + 1) - 1;

            // [start,end] is the range of outgoing labels in outlabels
            LL edge_l = outlabels->rank(start, c); // number of c-edges before the range
            LL edge_r = outlabels->rank(end+1, c); // number of c-edges up to the end of the range
            LL num_c_edges = edge_r - edge_l;

            if(num_c_edges == 0) return -1;

            l = indegs_rs.rank(indegs_ss0.select(C[kmer[i]] + edge_l + 1)) - 1; // -1: back to 0-indexing
            r = indegs_rs.rank(indegs_ss0.select(C[kmer[i]] + edge_r)) - 1; // -1: back to 0-indexing
        }

        assert(l == r); // Each k-mer should correspond to at most one node
        return l;        
    }

    // Returns the colex rank of the given kmer, or -1 if not found
    LL search(const string& kmer){
        assert(kmer.size() == k);
        return search(kmer.c_str());
    }
    
    // node ids go from 0 to (number of nodes) - 1
    LL walk(LL node, char c){

        if(node == -1) return -1;

        // Make sure the support structures point to the right vector (if the object is copied the addressed can change)
        set_supports(); 

        // Check if the node has an outedge with c. If no, return fail. If yes, we need to find the Wheeler rank 
        // of that outedge. Then we use the indegree vector to find the destination of that edge
        // See also:
        // Alanko, J., Gagie, T., Navarro, G., & Benkner, L. S. (2018). Tunneling on Wheeler Graphs.
        // arXiv preprint arXiv:1811.02457.
        
        LL start = outdegs_ss1(node+1) - node;

        LL end;
        if(node == n_nodes-1) end = outlabels->size()-1;
        else end = outdegs_ss1(node+2) - (node + 1) - 1;

        // The labels of outgoing nodes are in outlabels[start..end]
        if(outlabels->rank(end+1, c) - outlabels->rank(start, c) == 0){
            return -1; // Edge not found
        }

        LL edge_rank = C[c] + outlabels->rank(end+1,c); // \in [1, number of edges]
        LL pos_in_indegs = indegs_ss0(edge_rank);
        LL destination_rank = indegs_rs(pos_in_indegs); // \in [1, number of nodes]
        return destination_rank - 1; // -1: back to zero indexing of nodes
    }

    bool is_branching(LL node){ // todo: test
        // Make sure the support structures point to the right vector (if the object is copied the addressed can change)
        set_supports();
        
        LL i = outdegs_ss1(node+1); // Start of the unary number
        return i <= outdegs.size() - 3 && outdegs[i+1] == 0 && outdegs[i+2] == 0;
    }

    LL get_successor(LL node){
        if(node == -1) return -1;
        // Make sure the support structures point to the right vector (if the object is copied the addressed can change)
        set_supports();

        LL start = outdegs_ss1(node+1) - node;
        LL c = outlabels->char_at(start);
        return walk(node, c);
    }

    // There can be multiple incoming edges with the same label, but all of
    // them must have the same label. If there are multiple, returns an arbitrary choice of them.
    LL get_predecessor(LL node){
        if(node == -1) return -1;

        // Make sure the support structures point to the right vector (if the object is copied the addressed can change)
        set_supports(); 

        // Basically do what we do in the forward walk function, but in reverse
        
        LL edge_rank = indegs_ss1(node+1) - node; // zero-based rank

        // There are edge_rank edges before this edge. We want the largest c such
        // that C[c] <= edge_rank. This c is the incoming label to the node
        char c = 0;
        for(char a : alphabet){
            if(C[a] <= edge_rank) c = a;
        }
        LL rank_within_c_block = edge_rank - C[c]; // zero-based rank
        LL outedge_pos_in_outlabels = outlabels->select(rank_within_c_block+1, c);
        LL outedge_pos_in_outdegs = outdegs_ss0.select(outedge_pos_in_outlabels+1);
        return outdegs_rs.rank(outedge_pos_in_outdegs) - 1;
    }

    // Stores the outedge labels in the given vector and returns how many there are
    // DOES NOT RESIZE THE VECTOR. MAKE SURE IT HAS ENOUGH SPACE.
    // todo: unit test
    LL get_outlabels_from(LL node, vector<char>& answer){
        // Make sure the support structures point to the right vector (if the object is copied the addressed can change)
        set_supports();
        
        LL start = outdegs_ss1(node+1) - node;
        LL end;
        if(node == n_nodes-1) end = outlabels->size()-1;
        else end = outdegs_ss1(node+2) - (node + 1) - 1;

        LL count = 0;
        for(LL i = 0; i < end-start+1; i++){
            answer[i] = outlabels->char_at(start + i);
            count++;
        }
        return count;
    }

    // todo: unit test
    LL get_indegree(LL node){
        // Make sure the support structures point to the right vector (if the object is copied the addressed can change)
        set_supports();

        LL start = indegs_ss1(node+1) + 1;
        LL end;
        if(node == n_nodes-1) end = indegs.size()-1;
        else end = indegs_ss1(node+2) - 1;
        return end - start + 1;
    }

    // returns pair (node id, distance)
    pair<LL,LL> go_to_next_branching(LL node){
        if(is_branching(node)) return {node,0};
        if(get_number_of_edges() == get_number_of_nodes()){
            // The DBG is just a single cycle -> no branching nodes exist
            return {-1,-1};
        }

        // Make sure the support structures point to the right vector (if the object is copied the addressed can change)
        set_supports();

        // Now there exists at least one branching node. Since the input string is considered cyclic, all k-mers are reachable from all k-mers, se there exists a nearest branching node from the query node.

        // Reassigning the vector because the address of outdegs might have changed
        this->outdegs_ss1.set_vector(&(this->outdegs));
        
        LL start = outdegs_ss1(node+1) - node;
        pair<LL,LL> ans = go_to_next_branching(walk(node, outlabels->char_at(start)));
        return {ans.first, ans.second+1};
    }

    string to_string(){ // For debugging
        string bwt = get_outlabels();
        string bwt_spaced;

        LL next = 0;
        for(LL i = 0; i < outdegs.size(); i++){
            if(outdegs[i] == 1) bwt_spaced += " ";
            else{
                if(bwt[next] == BD_BWT_index<>::END) bwt_spaced += "*";
                else bwt_spaced += bwt[next];
                next++;
            }
        }
        
        stringstream ss;
        ss << bwt_spaced << "\n" << outdegs << "\n" << indegs << "\n";
        ss << "C = ";
        for(LL i = 0; i < C.size(); i++){
            if(i == 0 || C[i-1] != C[i]) ss<< C[i] << " "; 
        }
        return ss.str();
    }

    friend bool operator==(BOSS& boss1, BOSS& boss2);
};

bool operator==(BOSS& boss1, BOSS& boss2){
    // Returns true if the data is the same and the rank and
    // select structures give identical answers to all queries

    if(boss1.get_indegs() != boss2.get_indegs()) return false;
    if(boss1.get_outdegs() != boss2.get_outdegs()) return false;
    if(boss1.get_outlabels() != boss2.get_outlabels()) return false;
    if(boss1.get_C_array() != boss2.get_C_array()) return false;
    if(boss1.get_number_of_nodes() != boss2.get_number_of_nodes()) return false;
    if(boss1.get_k() != boss2.get_k()) return false;

    boss1.set_supports();
    boss2.set_supports();

    for(LL i = 0; i <= boss1.indegs.size(); i++){
        if(boss1.indegs_rs(i) != boss2.indegs_rs(i)) return false;
    }

    for(LL i = 1; i <= boss1.get_number_of_nodes(); i++){
        if(boss1.indegs_ss1(i) != boss2.indegs_ss1(i)) return false;
    }

    for(LL i = 1; i <= boss1.indegs.size() - boss1.get_number_of_nodes(); i++){
        if(boss1.indegs_ss0(i) != boss2.indegs_ss0(i)) return false;
    }

    for(LL i = 1; i <= boss1.get_number_of_nodes(); i++){
        if(boss1.outdegs_rs(i) != boss2.outdegs_rs(i)) return false;
    }

    for(LL i = 1; i <= boss1.get_number_of_nodes(); i++){
        if(boss1.outdegs_ss1(i) != boss2.outdegs_ss1(i)) return false;
    }

    for(LL i = 1; i <= boss1.get_number_of_nodes(); i++){
        if(boss1.outdegs_ss0(i) != boss2.outdegs_ss0(i)) return false;
    }

    if(boss1.alphabet != boss2.alphabet) return false;

    return true;

}

vector<LL> char_counts_to_C_array(vector<LL> counts){
    vector<LL> C(256); // Cumulative sum of counts

    // Compute cumulative sum of counts
    for(LL i = 0; i < C.size(); i++){
        C[i] = counts[i];
        if(i > 0) C[i] += C[i-1];
    }

    // Shift C to the right by one because that's how it's defined
    for(LL i = 256-1; i >= 0; i--){
        if(i == 0) C[i] = 0;
        else C[i] = C[i-1];
    }

    return C;

}

// Input: a file with all cyclic k-mers, one on each line, colex-sorted
// Splits the file by the last character of each kmer
vector<string> split_kmerfile(string kmerfile){
    map<char, ofstream*> out_streams;
    vector<string> filenames;

    string line;
    ifstream in(kmerfile);
    while(getline(in, line)){
        char last = line[line.size()-1];

        // If this is the first occurrence of character 'last'. create a new file for it
        if(out_streams.find(last) == out_streams.end()){
            string filename = get_temp_file_name("split");
            filenames.push_back(filename);
            ofstream* new_stream = new ofstream(filename);
            out_streams[last] = new_stream;
        }

        *(out_streams[last]) << line << "\n";
    }

    // Clean up
    for(auto keyval : out_streams){
        delete keyval.second;
    }

    return filenames;
}

// Input: a file with all cyclic (k+2)-mers, one on each line, colex-sorted
BOSS build_BOSS_from_kplus2_mers(string kmerfile, LL k){

    vector<bool> boss_indegrees; // Concatenation of unary representations. Indegree k -> 1 0^k
    vector<bool> boss_outdegrees;
    string boss_outlabels;
    vector<LL> boss_counts(256);

    // l-way merge. Put a pointer to the start of each file. 
    // Invariant: next largest element is on one of the pointed lines.
    // Store the pointed lines in a min-priority queue by the colex order of the middle k-mer.
    // Iteration:
    //   - Pop the minimum, push to the priority queue the next line from that file
    //   - While the minimum is still the same, pop that and push the next line
    //   - Record the indeg, outdeg and outlabels

    vector<string> kmerfiles = split_kmerfile(kmerfile);
    temp_file_manager.delete_file(kmerfile);

    vector<ifstream*> inputs;

    // Open files
    for(int64_t i = 0; i < kmerfiles.size(); i++) {
        ifstream* input = new ifstream(kmerfiles[i]);
        if(!input->good()){
            cerr << "Error opening file: " << kmerfiles[i] << endl;
            exit(1);
        }
        inputs.push_back(input);
    }

    auto cmp = [&](pair<string, int64_t> x, pair<string, int64_t> y) { 
        bool result = colex_compare(x.first.substr(1,k), y.first.substr(1,k));
        return result;
    };

    multiset<pair<string, int64_t>, decltype(cmp)> Q(cmp); // Priority queue: (kmer, file index).
    // Must be a multiset because a regular set will not add an element if it is equal
    // to another according to the comparison function.

    // Initialize priority queue
    string line;
    for(int64_t i = 0; i < inputs.size(); i++){
        getline(*(inputs[i]), line);
        assert(line.size() == k+2);
        Q.insert({line, i});
    }

    set<char> cur_outlabels;
    set<char> cur_inlabels;

    // Do the merge
    while(!Q.empty()){
        cur_outlabels.clear();
        cur_inlabels.clear();

        string kmer; LL stream_idx;
        std::tie(kmer, stream_idx) = *(Q.begin());
        Q.erase(Q.begin()); // pop

        cur_inlabels.insert(kmer[0]);
        cur_outlabels.insert(kmer.back());

        // Read next value from the file
        if(getline(*(inputs[stream_idx]), line)){
            assert(line.size() == k+2);
            Q.insert(make_pair(line, stream_idx));
        }

        while(!Q.empty() && Q.begin()->first.substr(1,k) == kmer.substr(1,k)){
            string next_kmer; LL next_stream_idx;
            std::tie(next_kmer, next_stream_idx) = *(Q.begin());
            cur_inlabels.insert(next_kmer[0]);
            cur_outlabels.insert(next_kmer.back());
            Q.erase(Q.begin()); // pop

            // Read next value from the file
            if(getline(*(inputs[next_stream_idx]), line)){
                assert(line.size() == k+2);
                Q.insert(make_pair(line, next_stream_idx));
            }
        }

        // Add the data to the BOSS arrays

        boss_outdegrees.push_back(1);
        for(char c : cur_outlabels){
            boss_outlabels += c;
            boss_outdegrees.push_back(0);
            boss_counts[c]++;
        }

        boss_indegrees.push_back(1);
        for(char c : cur_inlabels){
            (void)c; // Unused. Cast to void so that the compiler is happy
            boss_indegrees.push_back(0);
        }

    }

    for(ifstream* input : inputs) delete input; // Clean up
    for(string filename : kmerfiles) temp_file_manager.delete_file(filename); // Clean up

    // Put the data into the format that the BOSS constructor wants
    sdsl::bit_vector sdsl_indegs(boss_indegrees.size()), sdsl_outdegs(boss_outdegrees.size());
    for(LL i = 0; i < boss_indegrees.size(); i++) sdsl_indegs[i] = boss_indegrees[i];
    for(LL i = 0; i < boss_outdegrees.size(); i++) sdsl_outdegs[i] = boss_outdegrees[i];
    vector<LL> boss_C = char_counts_to_C_array(boss_counts);
    
    return BOSS(boss_outlabels, sdsl_indegs, sdsl_outdegs, boss_C, k);
}

BOSS build_BOSS_with_maps(const string& S, LL k){
    map<string, LL> k_counts; // k-mer counts. Keys are reverses of the k-mers so that the order is colexicographic
    map<string, LL> k_plus_1_counts; // k+1 mer counts. Keys are reverses of the k-mers so that the order is colexicographic

    set<char> alphabet;
    for(char c : S) alphabet.insert(c);

    for(LL i = 0; i < S.size(); i++){
        LL start = i;
        LL end = i + k - 1;

        {
            string rev_kmer;
            for(LL j = end; j >= start; j--) rev_kmer += S[j % S.size()];
            k_counts[rev_kmer]++;
        }
        
        end++;
        {
            string rev_kmer;
            for(LL j = end; j >= start; j--) rev_kmer += S[j % S.size()];
            k_plus_1_counts[rev_kmer]++;
        }
    }

    string outlabels;
    vector<bool> outdegrees; // Concatenation of unary representations. Indegree k -> 1 0^k
    vector<bool> indegrees;

    for(auto P : k_counts){
        // In colex order
        string rev_kmer = P.first;
        LL indegree = 0;
        LL outdegree = 0;
        for(char c : alphabet){
            if(k_plus_1_counts[c + rev_kmer] > 0){
                outlabels += c;
                outdegree++;
            }
            if(k_plus_1_counts[rev_kmer + c] > 0){
                indegree++;
            }
        }

        outdegrees.push_back(1);
        for(LL i = 0; i < outdegree; i++){
            outdegrees.push_back(0);
        }

        indegrees.push_back(1);
        for(LL i = 0; i < indegree; i++){
            indegrees.push_back(0);
        }
    }

    sdsl::bit_vector sdsl_indegs(indegrees.size()), sdsl_outdegs(outdegrees.size());
    for(LL i = 0; i < indegrees.size(); i++) sdsl_indegs[i] = indegrees[i];
    for(LL i = 0; i < outdegrees.size(); i++) sdsl_outdegs[i] = outdegrees[i];

    vector<LL> counts(256);
    for(char c : outlabels) counts[c]++;
    vector<LL> C = char_counts_to_C_array(counts);
    return BOSS(outlabels, sdsl_indegs, sdsl_outdegs, C, k);
}


sdsl::bit_vector to_sdsl(const vector<bool>& v){
    sdsl::bit_vector b(v.size());
    for(LL i = 0; i < v.size(); i++) b[i] = v[i];
    return b;
}

BOSS build_BOSS_with_bibwt(string S, LL k){
    write_log("Building the bibwt");
    std::reverse(S.begin(), S.end());
    // BWT and coBWT wlll be swapped. This is just because I have ready code
    // to mark kmers on the bwt but not on the coBWT. Now with these swapped I
    // can mark on the bwt which is actually the coBWT.
    BD_BWT_index<> bibwt((const uint8_t*)S.c_str());
    write_log("Marking k mers");
    sdsl::bit_vector k_marks = mark_kmers(bibwt,k);
    write_log("Marking k+1 mers");
    sdsl::bit_vector k_plus_1_marks = mark_kmers(bibwt,k+1);

    string outlabels;
    vector<bool> indegs;
    vector<bool> outdegs;
    set<char> distinct;
    vector<LL> counts(256);

    write_log("Building outlabels, outdegs, indegs and C");
    LL l = 0;
    LL r = 0;
    while(l < k_marks.size()){
        while(r < k_marks.size()-1 && k_marks[r+1] == 0) r++;
        // [l,r] is now the interval of a k-mer

        distinct.clear();
        LL indegree = 0;
        for(LL i = l; i <= r; i++){
            distinct.insert(bibwt.forward_bwt_at(i)); // It's actually the reverse bwt because we reversed the input
            if(k_plus_1_marks[i] == 1) indegree++;
        }

        outdegs.push_back(1);
        for(char c : distinct){
            outlabels += c;
            outdegs.push_back(0);
            counts[c]++;
        }

        indegs.push_back(1);
        for(LL i = 0; i < indegree; i++){
            indegs.push_back(0);
        }

        l = r+1;
        r = l;
    }

    BOSS boss(outlabels, to_sdsl(indegs), to_sdsl(outdegs), char_counts_to_C_array(counts), k);
    return boss;
}

BOSS build_BOSS_with_bibwt_from_fasta(string fastafile, LL k){
    string concat;
    FASTA_reader fr(fastafile);
    while(!fr.done()){
        Read_stream rs = fr.get_next_query_stream();
        concat += read_separator + rs.get_all();
    }
    return build_BOSS_with_bibwt(concat,k);
}


class Boss_Tester{
public:

    LL random_seed = 123674;

    class TestCase{
        public:
        string concat; // without BD_BWT_index<>::END
        vector<string> reads;
        set<string> kmers;
        set<string> k_plus1_mers;
        unordered_map<string,LL> kmer_to_node_id;
        set<char> alphabet;
        vector<string> colex_kmers;
        LL k;
    };

    // For KMC-based construction. KMC only supports the alphabet
    // {a,c,g,t,A,C,G,T}. It also turns all lower-case letters into upper case.
    class DNA_testcase{
        public:

        string concat;
        vector<string> reads;
        set<string> kmers;
        LL k;
    };

    Boss_Tester(){
        
    }

    void check_data_is_equal(BOSS& boss1, BOSS& boss2){
        assert(boss1 == boss2);
    }

    void test_serialization(TestCase& tcase){
        BOSS boss = build_BOSS_with_bibwt(tcase.concat, tcase.k); // Adds BD_BWT_index<>::END automatically
        boss.save_to_disk("test_out/");
        BOSS boss2;
        boss2.load_from_disk("test_out/");
        check_data_is_equal(boss, boss2);
    }

    vector<LL> get_neighbors(TestCase& tcase, LL id){
        string kmer = tcase.colex_kmers[id];
        vector<LL> neighbors;
        for(char c : tcase.alphabet){
            if(tcase.k_plus1_mers.count(kmer + c)){
                neighbors.push_back(tcase.kmer_to_node_id[kmer.substr(1)+c]);
            }
        }
        return neighbors;
    }

    vector<Boss_Tester::TestCase> generate_testcases(){
        srand(random_seed);
        vector<TestCase> testcases;
        LL n_reads = 10;
        for(LL read_length = 1; read_length <= 128; read_length *= 2){
            for(LL k = 1; k <= 256; k *= 2){
                k = min(k,(LL)254); // 254 + 2 = 256 is the biggest the EM sort can do
                TestCase tcase;
                vector<string> reads;
                string concat;
                for(LL i = 0; i < n_reads; i++) reads.push_back(get_random_dna_string(read_length,2));
                for(string S : reads) concat += read_separator + S;
                tcase.reads = reads;
                tcase.concat = concat;
                concat += BD_BWT_index<>::END;
                tcase.kmers = get_all_distinct_cyclic_kmers(concat,k);
                tcase.k_plus1_mers = get_all_distinct_cyclic_kmers(concat,k+1);
                tcase.colex_kmers = vector<string>(tcase.kmers.begin(), tcase.kmers.end());
                sort(tcase.colex_kmers.begin(), tcase.colex_kmers.end(), colex_compare);

                for(LL id = 0; id < tcase.colex_kmers.size(); id++)
                    tcase.kmer_to_node_id[tcase.colex_kmers[id]] = id;

                tcase.alphabet = get_alphabet(concat);
                tcase.k = k;
                testcases.push_back(tcase);
            }
        }
        return testcases;
    }

    void test_bibwt_construction(TestCase& tcase){

        BOSS boss1 = build_BOSS_with_bibwt(tcase.concat, tcase.k); // Adds BD_BWT_index<>::END automatically
        BOSS boss2 = build_BOSS_with_maps(tcase.concat + (char)(BD_BWT_index<>::END), tcase.k);
        check_data_is_equal(boss1, boss2);
        assert(tcase.kmers.size() == boss1.get_number_of_nodes());
        
    }

    void test_kplus2_construction(TestCase& tcase){

        // Write reads out in fasta
        string fastafile = temp_file_manager.get_temp_file_name("fasta");
        ofstream fasta_out(fastafile);
        
        for(string read : tcase.reads){
            fasta_out << ">\n" << read << "\n";
        }
        fasta_out.close();

        // List
        string kmers_outfile = temp_file_manager.get_temp_file_name("kmers_out");
        list_all_distinct_cyclic_kmers_in_external_memory(fastafile, kmers_outfile, tcase.k+2, 1, 2);

        // Sort
        string sorted_out = temp_file_manager.get_temp_file_name("kmers_sorted");

        auto colex_cmp = [](const char* x, const char* y){
            LL nx = strlen(x);
            LL ny = strlen(y);
            assert(nx != 0 && ny != 0);
            for(LL i = 0; i < min(nx,ny); i++){
                if(x[nx-1-i] < y[ny-1-i]) return true;
                if(x[nx-1-i] > y[ny-1-i]) return false;
            }
            // No mismatches -> the shorter string is smaller
            return nx < ny;
        };

        EM_sort(kmers_outfile, sorted_out, colex_cmp, 100, 4, 3, EM_LINES);

        BOSS boss1 = build_BOSS_from_kplus2_mers(sorted_out, tcase.k);
        BOSS boss2 = build_BOSS_with_maps(tcase.concat  + (char)(BD_BWT_index<>::END), tcase.k);
        cout << boss1.to_string() << endl << "==" << endl << boss2.to_string() << endl;
        check_data_is_equal(boss1, boss2);
        assert(tcase.kmers.size() == boss1.get_number_of_nodes());
        
    }

    void test_search(TestCase& tcase){
        // 1) If a k-mer is in concat, it must be found
        // 2) If a k-mer is not concat, it must not be found
        BOSS boss = build_BOSS_with_bibwt(tcase.concat, tcase.k); // Adds BD_BWT_index<>::END automatically

        // Check that all the k-mer in concat are found
        for(LL id = 0; id < tcase.colex_kmers.size(); id++){
            string kmer = tcase.colex_kmers[id]; // Must be found
            LL boss_id_searched = boss.search(kmer);
            assert(boss_id_searched != -1);
            assert(boss_id_searched == id);
        }

        // Check k-mers that are not found
        for(LL rep = 0; rep <= 1000; rep++){
            string kmer = get_random_string(tcase.k,20);
            if(tcase.kmers.count(kmer) == 0){
                // Not found
                assert(boss.search(kmer) == -1);
            }
        }
    }

    void test_walk(TestCase& tcase){
        BOSS boss = build_BOSS_with_bibwt(tcase.concat, tcase.k); // Adds BD_BWT_index<>::END automatically
        assert(boss.walk(-1,*tcase.alphabet.begin()) == -1); // Walk from -1 gives -1
        for(LL id = 0; id < tcase.colex_kmers.size(); id++){
            string kmer = tcase.colex_kmers[id];
            LL boss_id_searched = boss.search(kmer);
            assert(boss_id_searched == id);
            for(char c : tcase.alphabet){
                bool found = tcase.k_plus1_mers.count(kmer.substr(0) + c);
                LL boss_id = boss.walk(id, c);
                if(found){
                    assert(tcase.kmer_to_node_id[kmer.substr(1) + c] == boss_id);
                }
                else assert(boss_id == -1);
            }
        }
    }

    void test_get_successor(TestCase& tcase){
        BOSS boss = build_BOSS_with_bibwt(tcase.concat, tcase.k); // Adds BD_BWT_index<>::END automatically
        assert(boss.get_successor(-1) == -1); // Successor of -1 is - 1
        for(LL id = 0; id < tcase.colex_kmers.size(); id++){
            string kmer = tcase.colex_kmers[id];
            LL boss_id_searched = boss.search(kmer);
            assert(boss_id_searched == id);
            if(!boss.is_branching(id)){
                for(char c : tcase.alphabet){
                    LL succ_id = boss.walk(id, c);
                    if(succ_id != -1){
                        assert(succ_id == boss.get_successor(id));
                    }
                }
            }
        }
    }

    void test_get_predecessor(TestCase& tcase){
        BOSS boss = build_BOSS_with_bibwt(tcase.concat, tcase.k); // Adds BD_BWT_index<>::END automatically
        assert(boss.get_predecessor(-1) == -1); // Predecessor of -1 is -1
        for(LL id = 0; id < tcase.colex_kmers.size(); id++){
            string kmer = tcase.colex_kmers[id];
            LL boss_id_searched = boss.search(kmer);
            assert(boss_id_searched == id);
            for(char c : tcase.alphabet){
                LL succ_id = boss.walk(id, c);
                if(succ_id != -1 && boss.get_indegree(succ_id) == 1){
                    assert(id == boss.get_predecessor(succ_id));
                }
            }
        }
    }

    void test_is_branching(TestCase& tcase){
        BOSS boss = build_BOSS_with_bibwt(tcase.concat, tcase.k); // Adds BD_BWT_index<>::END automatically
        for(LL id = 0; id < boss.get_number_of_nodes(); id++){
            vector<LL> neighbors = get_neighbors(tcase,id);
            assert(neighbors.size() >= 1); // Since the input in considered cyclic, all k-mers should have at least one neighbor
            if(neighbors.size() == 1) assert(!boss.is_branching(id));
            else assert(boss.is_branching(id));
        }
    }

    void test_go_to_next_branching(TestCase& tcase){
        BOSS boss = build_BOSS_with_bibwt(tcase.concat, tcase.k); // Adds BD_BWT_index<>::END automatically
        for(LL id = 0; id < boss.get_number_of_nodes(); id++){
            if(tcase.kmers.size() == tcase.k_plus1_mers.size()){
                // The DBG is a single cycle -> no next branching node exists
                assert(boss.go_to_next_branching(id).first == -1);
            } else{
                LL node = id;
                LL distance = 0;
                while(true){
                    //cout << tcase.colex_kmers[node] << endl;
                    vector<LL> neighbors = get_neighbors(tcase,node);
                    assert(neighbors.size() >= 1); // Since the input in considered cyclic, all k-mers should have at least one neighbor
                    if(neighbors.size() >= 2) break;
                    else {
                        node = neighbors[0];
                        distance++;
                    }
                }
                assert(node == boss.go_to_next_branching(id).first);
                assert(distance == boss.go_to_next_branching(id).second);
            }
        }        
    }

};

void test_BOSS(){
    cerr << "Testing BOSS" << endl;
    disable_logging();

    Boss_Tester tester;
    for(Boss_Tester::TestCase tcase : tester.generate_testcases()){
        cerr << "Running testcase: " << "k = " << tcase.k << ", read_length = " << tcase.reads[0].size() << " n_reads = " << tcase.reads.size() << endl;
        cerr << "Testing (k+2)-mer construction" << endl;
        tester.test_kplus2_construction(tcase);
        cerr << "Testing bibwt construction" << endl;
        tester.test_bibwt_construction(tcase);
        cerr << "Testing serialization" << endl;
        tester.test_serialization(tcase);
        cerr << "Testing search" << endl;
        tester.test_search(tcase);
        cerr << "Testing walk" << endl;
        tester.test_walk(tcase);
        cerr << "Testing is_branching" << endl;
        tester.test_is_branching(tcase);
        cerr << "Testing go_to_next_branching" << endl;
        tester.test_go_to_next_branching(tcase);
        cerr << "Testing get_successor" << endl;
        tester.test_get_successor(tcase);
        cerr << "Testing get_predecessor" << endl;
        tester.test_get_predecessor(tcase);
    }

    enable_logging();
}
