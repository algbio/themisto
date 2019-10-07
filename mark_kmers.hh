#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include <cassert>
#include "test_tools.hh"
#include "BD_BWT_index.hh"

typedef int64_t LL;
using namespace std;

class Right_Maximal_Nodes_Iterator{
    // Iterates all right-maximal substrings including the empty string

    public:

    class Stack_frame{
    public:
        Interval_pair intervals;  // forward interval, reverse interval
        int64_t depth; // depth in the tree
        bool is_maxrep;
        Stack_frame(Interval_pair intervals, int64_t depth, bool is_maxrep) : intervals(intervals), depth(depth), is_maxrep(is_maxrep) {}
        Stack_frame(){}
    };
    
    std::stack<Stack_frame> iteration_stack;
    Stack_frame top;
    BD_BWT_index<>* index;
    typename BD_BWT_index<>::Interval_Data interval_data;
    
    Right_Maximal_Nodes_Iterator() {}
    Right_Maximal_Nodes_Iterator(BD_BWT_index<>* index) : index(index) {}
    
    void set_index(BD_BWT_index<>* index){ this->index = index; }
    
    
    Stack_frame get_top(){ return top; }
    
    
    void init(){
        // Make space
        interval_data.symbols.resize(index->get_alphabet().size());
        interval_data.ranks_start.resize(index->get_alphabet().size());
        interval_data.ranks_end.resize(index->get_alphabet().size());
        
        // Clear the stack
        while(!iteration_stack.empty()) iteration_stack.pop();
        
        // Push the empty string 
        iteration_stack.push(Stack_frame(Interval_pair(0,index->size()-1,0,index->size()-1),0,true));

    }
    
    // After calling this, the member variable 'top' contains the stack frame of the next string.
    bool next(){

         if(iteration_stack.empty()) return false;
         
         top = iteration_stack.top(); iteration_stack.pop();
         
         index->compute_bwt_interval_data(top.intervals.forward, interval_data);
         
         // Iterate alphabet in reverse lexicographic order, so the smallest is pushed to the
         // stack the last, so the iteration is done in lexicographic order.
         // Push all right-maximal children.
         for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
            Interval_pair I2 = index->left_extend(top.intervals,interval_data,i);
            if(I2.forward.size() != 0 && index->is_right_maximal(I2)){
                iteration_stack.push(Stack_frame(I2, top.depth+1, index->is_right_maximal(I2)));
            }
        }
        
        return true;
    }
};

sdsl::bit_vector mark_kmers(BD_BWT_index<>& bibwt, LL k){
    // Algorithm: we want to mark the first leaf of every ST node whole depth is at least k
    // but the depth of the parent is less than k. So: iterate all internal ST nodes of depth
    // less than k and at each, mark the first leaf of every child.

    Right_Maximal_Nodes_Iterator it(&bibwt);
    sdsl::bit_vector marks(bibwt.size(),0);
    Interval_pair child; // Allocate this here so we can reuse to avoid slow heap access inside the loop
    it.init();
    while(it.next()){
        if(it.top.depth < k){
            bibwt.compute_rev_bwt_interval_data(it.top.intervals.reverse, it.interval_data);
            for(int64_t i = it.interval_data.n_distinct_symbols-1; i >= 0; i--){
                child = bibwt.right_extend(it.top.intervals,it.interval_data,i);
                marks[child.forward.left] = 1;
            }
        }
    }
    return marks;
}

sdsl::bit_vector mark_kmers_brute(string S, LL k){
    vector<string> suffixes;
    for(LL i = 0; i < S.size(); i++){
        suffixes.push_back(S.substr(i));
    }
    sort(suffixes.begin(), suffixes.end());

    sdsl::bit_vector marks(S.size(),0);
    marks[0] = 1;
    for(LL i = 1; i < S.size(); i++){
        if(suffixes[i].substr(0,k) != suffixes[i-1].substr(0,k)){
            marks[i] = 1;
        }
    }
    return marks;
}


void mark_kmers_tests(){
    cerr << "Running mark kmers tests" << endl;

    // todo: large test case

    //string S = "bbabbbbabbbbababbaabaaaabbabaa#abbabaaaababaaaabbabaabaabbbaa#bababaabbbababbbbaaabbabbbabaa#bbbababbbaabbbaaaaaaabbaabbaab#bbaabbaabbabababbabbbbbabbbbab#ababbaabababbaabbaaabbbabbabaa#aabbbababbbbaabbabbababbbaabba#baaabaaaabbbaaabbaaaaaabbbbbab#bbbbbabbababbabbaabaababababab#ababbababbaabbaaaabbbbabaabbbb#aababaabaaaaabbaabbbaaaabbaaaa#aaababbbabbbababbababbaaaabaab#aababbaaaaabbaaabaaaaaabaabbab#ababbbababaaababbbbaabaaaaabba#bbbbabaaaabaaaaabbaababbaababa#aababbbaabababbababbabbaabaabb#aaaaabbbbbaababaaabaabbbbbbbba#baabaaabbabaaaababbbbbbababaab#abbbabbaababbabbbaabbbaabbbbbb#aaabbaaaaabbaabbabbaabaabbaaab#";
    //BD_BWT_index<> bibwt((const uint8_t*)(S.c_str()));
    //cout << mark_kmers(bibwt,2) << endl;
    //cout << mark_kmers_brute(S,2) << endl;
    //exit(0);

    LL maxlen = 10;
    for(string S : all_binary_strings_up_to(maxlen)){
        for(LL k = 1; k < 4; k++){
            BD_BWT_index<> bibwt((const uint8_t*)(S.c_str()));
            assert(mark_kmers(bibwt, k) == mark_kmers_brute(S + '\x01',k));
        }
    }
    cerr << "mark kmers for all binary string up to length " << maxlen << " ok" << endl;

    for(LL rep = 0; rep < 100; rep++){
        string S = get_random_dna_string(40, 3);
        cout << S << endl;
        LL k = rand() % 5 + 1;
        BD_BWT_index<> bibwt((const uint8_t*)(S.c_str()));
        assert(mark_kmers(bibwt, k) == mark_kmers_brute(S + '\x01',k)); // Add the bibwt end sentinel to brute
    }
}