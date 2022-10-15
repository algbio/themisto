#include <string>
#include <sstream>
#include "pseudoalign.hh"
#include "WorkDispatcher.hh"

using namespace std;
using namespace sbwt;


template <typename T>
string vec_to_string(const vector<T>& v){
    stringstream ss;
    for(T x : v) ss << x << " ";
    return ss.str();
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

int64_t fast_int_to_string(int64_t x, char* buffer){
    // Fast manual integer-to-string conversion

    LL i = 0;
    // Write the digits in reverse order (reversed back at the end)
    if(x == -1){
        buffer[0] = '1';
        buffer[1] = '-';
        i = 2;
    } else if(x == 0){
        buffer[0] = '0';
        i = 1;
    } else{
        while(x > 0){
            buffer[i++] = '0' + (x % 10);
            x /= 10;
        }
    }
    std::reverse(buffer, buffer + i);
    buffer[i] = '\0';
    return i;
}