#include <string>
#include "pseudoalign.hh"

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

void pseudoalign(const plain_matrix_sbwt_t& SBWT, const Coloring& coloring, int64_t n_threads, std::string inputfile, std::string outputfile, bool reverse_complements, int64_t buffer_size, bool gzipped, bool sorted_output){

    SeqIO::Reader<> reader(inputfile);
    LL string_id = 0;
    while(true) { 
        LL len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        vector<int64_t> colex_ranks = SBWT.streaming_search(reader.read_buf, len);
        Color_Set intersection;
        bool first_nonempty_found = false;
        for(LL colex : colex_ranks){
            if(colex >= 0){ // k-mer found Found
                Color_Set cs = coloring.get_color_set(colex);
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

        // Todo: write to output stream with buffering and parallel safety etc
        cout << string_id << " ";
        for(color_t x : intersection.get_colors_as_vector()){
            cout << x <<  " ";
        }
        cout << "\n";
        string_id++;
    }

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