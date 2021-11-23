#pragma once

#include "../globals.hh"
#include "../tests/setup_tests.hh"
#include <gtest/gtest.h>
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

TEST(KMER, basic){
    for(LL len = 1; len <= 255; len++){
        string S = debug_test_get_random_DNA_string(len);
        Kmer<255> kmer(S);

        // Check that the constructor worked right
        ASSERT_TRUE(kmer.get_k() == S.size());
        for(LL i = 0; i < len; i++){
            ASSERT_TRUE(kmer.get(i) == S[i]);
        }

        // Test copy
        Kmer<255> kmer_copy = kmer.copy();
        ASSERT_TRUE(kmer_copy == kmer);

        // Edit characters randomly
        for(LL i = 0; i < 1000; i++){
            LL idx = rand() % S.size();
            char c = get_random_DNA_char();
            kmer.set(idx, c);
            S[idx] = c;
            ASSERT_TRUE(kmer.get(idx) == S[idx]);
        }

        // Check that the edits worked out right
        for(LL i = 0; i < len; i++){
            ASSERT_TRUE(kmer.get(i) == S[i]);
        }

        // Test drop left
        Kmer<255> left = kmer.copy().dropleft();
        ASSERT_TRUE(left.get_k() == len-1);
        for(LL i = 0; i < len-1; i++){
            ASSERT_TRUE(left.get(i) == S[i+1]);
        }

        // Test append left
        Kmer<255> left_append = left.copy().appendleft(kmer.first());
        ASSERT_TRUE(left_append.get_k() == len);
        ASSERT_TRUE(left_append == kmer);

        // Test drop right
        Kmer<255> right = kmer.copy().dropright();
        ASSERT_TRUE(right.get_k() == len-1);
        for(LL i = 0; i < len-1; i++){
            ASSERT_TRUE(right.get(i) == S[i]);
        }

        // Test append right
        Kmer<255> right_append = right.copy().appendright(kmer.last());
        ASSERT_TRUE(right_append.get_k() == len);
        ASSERT_TRUE(right_append == kmer);
    }
}

TEST(KMER, serialization){
    const string S = debug_test_get_random_DNA_string(255);
    const Kmer<255> kmer(S);

    // Test serialization to file
    string filename = get_temp_file_manager().create_filename("kmer-serialization");
    ofstream out(filename, ios_base::binary);
    kmer.serialize(out);
    out.flush();
    ifstream in(filename, ios_base::binary);
    Kmer<255> loaded;
    loaded.load(in);
    ASSERT_TRUE(kmer == loaded);
    ASSERT_TRUE(loaded.to_string() == S);
    ASSERT_TRUE(loaded.get_k() == S.size());

    // Test serialization to a char buffer
    char buffer[Kmer<255>::size_in_bytes()];
    kmer.serialize(buffer);
    Kmer<255> loaded2;
    loaded2.load(buffer);
    ASSERT_TRUE(kmer == loaded2);
    ASSERT_TRUE(loaded2.to_string() == S);
    ASSERT_TRUE(loaded2.get_k() == S.size());    
}

TEST(KMER, colex){
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
                ASSERT_TRUE(kmers[i] == kmers[j]);
                ASSERT_TRUE((kmers[i] < kmers[j]) == false);
            } else{
                ASSERT_TRUE(colex_compare(strings[i], strings[j]) == (kmers[i] < kmers[j]));
            }
        }
    }
}