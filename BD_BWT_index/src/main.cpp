#include <iostream>
#include <string>
#include "BD_BWT_index.hh"
#include "Iterators.hh"
#include <cassert>
#include <set>

using namespace std;

static int64_t modulo(int64_t x, int64_t m){
    int64_t r = x % m;
    if(r < 0) r += m;
    return r;
}

set<string> get_right_maximal_substrings(const string& s){
    int n = s.size();
    set<string> result;
    set<char> alphabet(s.begin(),s.end());
    set<string> substrings;
    substrings.insert("");
    for(int k = 1; k <= n; k++){
        for(int i = 0; i <= n - k; i++){
            substrings.insert(s.substr(i,k));
        }
    }
    for(const string& x : substrings){
        int nExtensions = 0; // Number of ways to extend x right. End of string counts a distinct extension.
        for(char c : alphabet){
            if(substrings.count(x + c) == 1)
                nExtensions++;
        }
        if(s.substr(n-x.size()) == x) nExtensions++; // x is at the right end of s
        // Note subtr: "If this is equal to the string length, the function returns an empty string."
        
        if(nExtensions >= 2)
            result.insert(x);
    }
    return result;
}

vector<string> all_binary_strings_up_to(int k){
    vector<string> ans;
    ans.push_back("");
    for(int length = 1; length <= k; length++){
        for(int mask = 0; mask < (1 << length); mask++){
            string s = "";
            for(int i = 0; i < length; i++){
                if(mask & (1 << i)) s += 'a';
                else s += 'b';
            }
            ans.push_back(s);
        }
    }
    return ans;
}


bool test_suffix_link_tree_iteration(BD_BWT_index<sdsl::bit_vector>& index, string& s){
    BD_BWT_index_iterator<sdsl::bit_vector> it(&index);
    set<string> labels;
    while(it.next()){
        string x(it.label.rbegin(), it.label.rend());
        labels.insert(x);
    }
    
    return (labels == get_right_maximal_substrings(s));       
}

bool test_backward_step(BD_BWT_index<sdsl::bit_vector>& index, const string& s){
    BD_BWT_index_iterator<sdsl::bit_vector> it(&index);
    string s_with_end = s + (char)BD_BWT_index<sdsl::bit_vector>::END;
    int64_t n = s_with_end.size();
    int64_t lex_rank = 0; // Lex rank of the "empty suffix"
    for(int64_t i = 0; i < n*3; i++){
        int64_t text_pos = modulo(n-1-i, n);
        if(index.forward_bwt_at(lex_rank) != s_with_end[modulo(text_pos-1,n)]) return false;
        lex_rank = index.backward_step(lex_rank);
    }
    return true;
}

bool test_forward_step(BD_BWT_index<sdsl::bit_vector>& index, const string& s){
    BD_BWT_index_iterator<sdsl::bit_vector> it(&index);
    string s_with_end = s + (char)BD_BWT_index<sdsl::bit_vector>::END;
    int64_t n = s_with_end.size();
    int64_t colex_rank = 0; // Colexicographic rank of the "empty prefix"
    for(int64_t i = 0; i < n*3; i++){
        int64_t text_pos = modulo(n-1+i, n);
        if(index.backward_bwt_at(colex_rank) != s_with_end[modulo(text_pos+1,n)]) return false;
        colex_rank = index.forward_step(colex_rank);
    }
    return true;
}

int main(int argc, char** argv){
    
    vector<string> test_set = all_binary_strings_up_to(10);
    for(auto& s : test_set){
        if(s == "") continue;
        cout << s << endl;
        BD_BWT_index<> index((const uint8_t*)s.c_str());
        assert(test_suffix_link_tree_iteration(index,s));
        assert(test_backward_step(index,s));
        assert(test_forward_step(index,s));
    }
    
    cerr << "All tests OK" << endl;
    
}