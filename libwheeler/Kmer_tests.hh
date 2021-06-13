#pragma once

#include "../globals.hh"
#include "Kmer.hh"

char get_random_DNA_char(){
    LL r = rand() % 4;
    if(r == 0) return 'A';
    else if(r == 1) return 'C';
    else if(r == 2) return 'G';
    else return 'T';
}

string debug_test_get_random_DNA_string(LL len){
    string S;
    for(LL i = 0; i < len; i++){
        S += get_random_DNA_char();
    }
    return S;
}

void Kmer_basic_tests(){
    for(LL len = 1; len <= 255; len++){
        string S = debug_test_get_random_DNA_string(len);
        Kmer<255> kmer(S);

        // Check that the constructor worked right
        assert(kmer.get_k() == S.size());
        for(LL i = 0; i < len; i++){
            assert(kmer.get(i) == S[i]);
        }

        // Test copy
        Kmer<255> kmer_copy = kmer.copy();
        assert(kmer_copy == kmer);

        // Edit characters randomly
        for(LL i = 0; i < 1000; i++){
            LL idx = rand() % S.size();
            char c = get_random_DNA_char();
            kmer.set(idx, c);
            S[idx] = c;
            assert(kmer.get(idx) == S[idx]);
        }

        // Check that the edits worked out right
        for(LL i = 0; i < len; i++){
            assert(kmer.get(i) == S[i]);
        }

        // Test drop left
        Kmer<255> left = kmer.copy().dropleft();
        assert(left.get_k() == len-1);
        for(LL i = 0; i < len-1; i++){
            assert(left.get(i) == S[i+1]);
        }

        // Test append left
        Kmer<255> left_append = left.copy().appendleft(kmer.first());
        assert(left_append.get_k() == len);
        assert(left_append == kmer);

        // Test drop right
        Kmer<255> right = kmer.copy().dropright();
        assert(right.get_k() == len-1);
        for(LL i = 0; i < len-1; i++){
            assert(right.get(i) == S[i]);
        }

        // Test append right
        Kmer<255> right_append = right.copy().appendright(kmer.last());
        assert(right_append.get_k() == len);
        assert(right_append == kmer);
    }
    cerr << "Kmer copy, access, modify, dropleft, appendleft, dropright and appendright tests passed" << endl;
}

void serialization_test(){
    string S = debug_test_get_random_DNA_string(255);
    Kmer<255> kmer(S);
    ofstream out("serialization_test.bin", ios_base::binary);
    kmer.serialize(out);
    out.flush();
    ifstream in("serialization_test.bin", ios_base::binary);
    Kmer<255> loaded;
    loaded.load(in);
    assert(kmer == loaded);
    assert(kmer.to_string() == S);
    assert(kmer.get_k() == S.size());
    cerr << "Kmer serialization test passed" << endl;
}

void Kmer_colex_tests(){
    vector<string> strings;
    vector<Kmer<255>> kmers;

    // Generate random k-mers
    for(int i = 0; i < 20; i++){
        strings.push_back(debug_test_get_random_DNA_string(rand() % 256)); // Random length, alphabet size 4
    }

    // Add empty string
    strings.push_back("");
    strings.push_back(""); // Another so we have empty vs empty comparison

    // Some max-length strings for good measure
    for(LL i = 0; i < 10; i++){
        strings.push_back(debug_test_get_random_DNA_string(255));
    }

    // Add some strings that have a long shared suffix
    for(LL i = 0; i < 40; i++){
        string suffix = debug_test_get_random_DNA_string(100 + i);
        strings.push_back(suffix); // Exact suffix of another
        strings.push_back(suffix); // Lets put in another so we have a long exact match
        strings.push_back("AAAAAAAAAAAAAA" + suffix); // Prepend some small mismatches
        strings.push_back("TTTTTTTTTTTTTT" + suffix); // Prepend some large mismatches
    }

    // Build the k-mer objects
    for(string S : strings) kmers.push_back(Kmer<255>(S));

    // Check all comparisons
    for(int i = 0; i < strings.size(); i++){
        for(int j = 0; j < strings.size(); j++){
            if(strings[i] == strings[j]){
                assert(kmers[i] == kmers[j]);
                assert((kmers[i] < kmers[j]) == false);
            } else{
                assert(colex_compare(strings[i], strings[j]) == (kmers[i] < kmers[j]));
            }
        }
    }

    cerr << "Kmer colex comparison test passed" << endl;
}

void kmer_tests(){
    Kmer_basic_tests();
    Kmer_colex_tests();
    serialization_test();
}