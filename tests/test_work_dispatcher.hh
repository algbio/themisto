#pragma once

#include <gtest/gtest.h>
#include <vector>
#include <unordered_map>
#include "setup_tests.hh"
#include "globals.hh"
#include "WorkDispatcher.hh"
#include "sbwt/SeqIO.hh"

class DispatcherConsumerTestCallback : DispatcherConsumerCallback{

public:

    bool finished_flag = false;
    vector<string> received_strings;
    vector<int64_t> received_string_ids;
    LL busywork = 0;

    virtual void callback(const char* S, LL S_size, int64_t string_id){
        received_strings.push_back(string(S,S_size));
        received_string_ids.push_back(string_id);
        for(LL i = 0; i < 100000; i++) busywork++; // Do some "work"
    }

    virtual void finish(){
        finished_flag = true;
    }

    virtual ~DispatcherConsumerTestCallback() {} 
};

TEST(WORK_DISPATCHER, basic_test){
    // void run_dispatcher(vector<DispatcherConsumerCallback*>& callbacks, Sequence_Reader_Buffered& sr, LL buffer_size);
    vector<string> seqs;

    // Create 10MB of sequence data
    for(LL i = 0; i < 1e5; i++){
        seqs.push_back(get_random_dna_string(100,4));
    }

    string fastafile = get_temp_file_manager().create_filename("",".fna");
    write_as_fasta(seqs, fastafile);

    sbwt::SeqIO::Reader<> sr(fastafile);
    LL buffer_size = 1024; // Small buffer so that not just one thread gets all the work
    
    vector<DispatcherConsumerCallback*> callbacks;
    vector<DispatcherConsumerTestCallback*> callbacks_no_cast;
    for(LL i = 0; i < 4; i++){
        DispatcherConsumerTestCallback* cb = new DispatcherConsumerTestCallback();
        callbacks.push_back((DispatcherConsumerCallback*)cb);
        callbacks_no_cast.push_back(cb);
    }

    run_dispatcher(callbacks, sr, buffer_size);

    vector<pair<LL, string> > pairs; // pairs (string id, string)
    for(LL i = 0; i < 4; i++){
        DispatcherConsumerTestCallback* cb = callbacks_no_cast[i];
        ASSERT_EQ(cb->received_string_ids.size(), cb->received_strings.size());
        for(LL j = 0; j < cb->received_strings.size(); j++){
            pairs.push_back({cb->received_string_ids[j], cb->received_strings[j]});
        }
    }

    ASSERT_EQ(pairs.size(), seqs.size());
    std::sort(pairs.begin(), pairs.end());
    for(LL i = 0; i < pairs.size(); i++){
        ASSERT_EQ(pairs[i].first, i);
        ASSERT_EQ(pairs[i].second, seqs[i]);
    }

    for(LL i = 0; i < 4; i++)
        delete callbacks[i];

}
