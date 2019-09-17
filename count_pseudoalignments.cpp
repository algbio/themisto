#include <string>
#include <unordered_map>
#include <iostream>
#include <vector>
#include "globals.hh"

using namespace std;

int main(){
    // Reads input in the output format of the pseudoalign program from stdin, i.e.
    // On each line there is a list of numbers. The first is the sequence
    // id of the query read, and the rest are sequence ids of the references
    // that this query pseudoaligns to.

    // Outputs lines with pairs (reference id, number of reads pseudoaligned to this reference)

    unordered_map<LL, LL> counts;
    string line;
    while(getline(cin, line)){
        vector<LL> numbers = parse_tokens<LL>(line);
        for(LL i = 1; i < numbers.size(); i++){ // Ignore first token becaue it's the read id
            counts[numbers[i]]++;
        }
    }

    for(auto keyvalue : counts){
        LL ref_id = keyvalue.first;
        LL count = keyvalue.second;
        cout << ref_id << " " << count << endl;
    }

}