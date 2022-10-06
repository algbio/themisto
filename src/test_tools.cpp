#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include <cassert>
#include <set>
#include <unordered_map>
#include <map>
#include "globals.hh"
#include "stdlib_printing.hh"
#include "test_tools.hh"
#include <cassert>

using namespace std;
typedef int64_t LL;

set<char> get_alphabet(string S){
    set<char> ans;
    for(char c: S) ans.insert(c);
    return ans;
}

set<string> get_all_distinct_kmers(string S, int64_t k){
    set<string> kmers;
    for(LL i = 0; i < S.size()-k+1; i++){
        kmers.insert(S.substr(i,k));
    }
    return kmers;
}

set<string> get_all_distinct_cyclic_kmers(string S, int64_t k){
    set<string> kmers;
    for(int64_t i = 0; i < S.size(); i++){
        string kmer;
        for(int64_t j = 0; j < k; j++){
            kmer += S[(i+j) % S.size()];
        }
        kmers.insert(kmer);
    }
    return kmers;
}


set<string> get_all_distinct_cyclic_kmers(vector<string>& A, int64_t k){
    throw std::runtime_error("ERROR: OLD CODE");
    /*
    string concat;
    for(string read : A){
        concat += read_separator + read;
    }
    concat += '\x01'; // bibwt end sentinel

    return get_all_distinct_cyclic_kmers(concat,k);
    */
}

vector<string> get_all_kmers(string& S, int64_t k){
    vector<string> kmers;
    for(int64_t i = 0; i < S.size()-k+1; i++){
        kmers.push_back(S.substr(i,k));
    }
    return kmers;
}

vector<string> all_binary_strings_up_to(int64_t k){ // For testing
    vector<string> ans;
    for(int64_t length = 1; length <= k; length++){
        for(int64_t mask = 0; mask < (1 << length); mask++){
            string s = "";
            for(int64_t i = 0; i < length; i++){
                if(mask & (1 << i)) s += 'A';
                else s += 'C';
            }
            ans.push_back(s);
        }
    }
    return ans;
}

string get_random_dna_string(int64_t length, int64_t alphabet_size){ // For testing
    string s;
    assert(alphabet_size >= 1 && alphabet_size <= 4);
    char alphabet[4] = {'A','T','G','C'};
    for(int64_t i = 0; i < length; i++){
        s.push_back(alphabet[rand() % alphabet_size]);
    }
    return s;
}

string get_random_string(int64_t length, int64_t alphabet_size){ // For testing
    string s;
    for(int64_t i = 0; i < length; i++){
        int64_t r = rand() % alphabet_size;
        s += 'a' + r;
    }
    return s;
}

vector<string> get_sorted_suffixes(string S){
    vector<string> suffixes;
    for(int64_t i = 0; i < S.size(); i++){
        suffixes.push_back(S.substr(i));
    }
    sort(suffixes.begin(), suffixes.end());
    return suffixes;
}

void write_as_fasta(vector<string>& seqs, string fasta_filename){
    sbwt::throwing_ofstream out(fasta_filename);
    for(string& S : seqs) out << ">\n" << S << "\n";
}

void write_as_fastq(vector<string>& seqs, string fastq_filename){
    sbwt::throwing_ofstream out(fastq_filename);
    for(string& S : seqs){
        out << "@\n" << S << "\n+\n" << string(S.size(), 'I') << "\n";
    }
}

static char incoming_label(const sbwt::plain_matrix_sbwt_t& SBWT, int64_t node){
    if(node < SBWT.get_C_array()[0]) return '$';
    else if(node < SBWT.get_C_array()[1]) return 'A';
    else if(node < SBWT.get_C_array()[2]) return 'C';
    else if(node < SBWT.get_C_array()[3]) return 'G';
    else return 'T';
}


// Include dummy nodes
vector<string> dump_node_labels(sbwt::plain_matrix_sbwt_t& SBWT){
    vector<sdsl::select_support_mcl<>> select_supports(4);
    sdsl::util::init_support(select_supports[0], &SBWT.get_subset_rank_structure().A_bits);
    sdsl::util::init_support(select_supports[1], &SBWT.get_subset_rank_structure().C_bits);
    sdsl::util::init_support(select_supports[2], &SBWT.get_subset_rank_structure().G_bits);
    sdsl::util::init_support(select_supports[3], &SBWT.get_subset_rank_structure().T_bits);

    vector<string> labels;
    LL k = SBWT.get_k();
    for(LL i = 0; i < SBWT.number_of_subsets(); i++){
        string label;
        LL node = i;
        for(LL j = 0; j < k; j++){
            char c = incoming_label(SBWT, node);
            label.push_back(c);
            if(c == '$') 
                continue;
            if(c == 'A')
                node = select_supports[0].select(node - SBWT.get_C_array()[0] + 1);
            if(c == 'C')
                node = select_supports[1].select(node - SBWT.get_C_array()[1] + 1);
            if(c == 'G')
                node = select_supports[2].select(node - SBWT.get_C_array()[2] + 1);
            if(c == 'T')
                node = select_supports[3].select(node - SBWT.get_C_array()[3] + 1);
        }

        std::reverse(label.begin(), label.end());
        labels.push_back(label);
    }
    return labels;
}