#include <string>
#include "pseudoalign.hh"
#include <sstream>

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

void pseudoalign(const plain_matrix_sbwt_t& SBWT, const Coloring& coloring, int64_t n_threads, std::string inputfile, std::string outputfile, bool reverse_complements, int64_t buffer_size, bool gzipped, bool sorted_output){

    SeqIO::Reader<> reader(inputfile);
    LL string_id = 0;
    sbwt::throwing_ofstream* out = nullptr;
    if(outputfile != "") out = new sbwt::throwing_ofstream(outputfile);
    while(true) { 
        LL len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        vector<int64_t> colex_ranks = SBWT.streaming_search(reader.read_buf, len);
        vector<int64_t> rc_colex_ranks;
        if(reverse_complements){
            reverse_complement_c_string(reader.read_buf, len);
            rc_colex_ranks = SBWT.streaming_search(reader.read_buf, len);
            reverse_complement_c_string(reader.read_buf, len); // Undo reverse complement
        }
        LL n_kmers = colex_ranks.size();
        Color_Set intersection;
        bool first_nonempty_found = false;
        for(LL i = 0; i < n_kmers; i++){
            if(colex_ranks[i] >= 0 || (reverse_complements && rc_colex_ranks[n_kmers-1-i] >= 0)){ // k-mer found
                Color_Set cs; // Empty
                if(colex_ranks[i] >= 0) cs = coloring.get_color_set(colex_ranks[i]);
                cout << "cs " << vec_to_string(cs.get_colors_as_vector()) << endl;
                if(reverse_complements && rc_colex_ranks[n_kmers-1-i] >= 0){
                    Color_Set cs_rc = coloring.get_color_set(rc_colex_ranks[n_kmers-1-i]);
                    cout << "cs_rc " << vec_to_string(cs_rc.get_colors_as_vector()) << endl;
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

        stringstream output_buffer;
        output_buffer << string_id << " ";
        for(color_t x : intersection.get_colors_as_vector()){
            output_buffer << x <<  " ";
        }
        output_buffer << "\n";

        if(out == nullptr) cout << output_buffer.str();
        else out->stream << output_buffer.str();

        string_id++;

        // Todo: write to output stream with buffering and parallel safety etc
    }

    delete out;

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