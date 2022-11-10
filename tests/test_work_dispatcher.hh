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
    vector<int64_t> received_metadata;
    int64_t busywork = 0;

    virtual void callback(const char* S, int64_t S_size, int64_t string_id, std::array<uint8_t, 8> metadata){
        received_strings.push_back(string(S,S_size));
        received_string_ids.push_back(string_id);
        received_metadata.push_back(*reinterpret_cast<int64_t*>(metadata.data()));
        for(int64_t i = 0; i < 100000; i++) busywork++; // Do some "work"
    }

    virtual void finish(){
        finished_flag = true;
    }

    virtual ~DispatcherConsumerTestCallback() {} 
};

// Color stream from an in-memory vector
class In_Memory_Int_Stream : public Metadata_Stream{

private:

    vector<int64_t> ints;
    int64_t idx = 0; // Current index
    std::array<uint8_t, 8> dummy;

public:

    In_Memory_Int_Stream(vector<int64_t> ints) : ints(ints) {}

    virtual std::array<uint8_t, 8> next(){
        if(idx == ints.size()) return dummy; // Done

        std::array<uint8_t, 8> ret;
        int64_t* ptr = (int64_t*)ret.data(); // Interpret as int64_t
        *ptr = ints[idx++];
        
        return ret;
    }

};

TEST(WORK_DISPATCHER, basic_test){
    vector<string> seqs;
    vector<int64_t> metadata;

    // Create 10MB of sequence data
    for(int64_t i = 0; i < 1e5; i++){
        seqs.push_back(get_random_dna_string(100,4));
        metadata.push_back(100 + i);
    }

    string fastafile = get_temp_file_manager().create_filename("",".fna");
    write_as_fasta(seqs, fastafile);

    sbwt::SeqIO::Reader<> sr(fastafile);
    int64_t buffer_size = 1024; // Small buffer so that not just one thread gets all the work
    
    vector<DispatcherConsumerCallback*> callbacks;
    vector<DispatcherConsumerTestCallback*> callbacks_no_cast;
    for(int64_t i = 0; i < 4; i++){
        DispatcherConsumerTestCallback* cb = new DispatcherConsumerTestCallback();
        callbacks.push_back((DispatcherConsumerCallback*)cb);
        callbacks_no_cast.push_back(cb);
    }

    In_Memory_Int_Stream metadata_stream(metadata);
    run_dispatcher(callbacks, sr, &metadata_stream, buffer_size);

    vector<tuple<int64_t, string, int64_t> > received; // pairs (string id, string, metadata)
    for(int64_t i = 0; i < 4; i++){
        DispatcherConsumerTestCallback* cb = callbacks_no_cast[i];
        ASSERT_EQ(cb->received_string_ids.size(), cb->received_strings.size());
        for(int64_t j = 0; j < cb->received_strings.size(); j++){
            received.push_back({cb->received_string_ids[j], cb->received_strings[j], cb->received_metadata[j]});
        }
    }

    ASSERT_EQ(received.size(), seqs.size());
    std::sort(received.begin(), received.end());
    for(int64_t i = 0; i < received.size(); i++){
        ASSERT_EQ(get<0>(received[i]), i);
        ASSERT_EQ(get<1>(received[i]), seqs[i]);
        ASSERT_EQ(get<2>(received[i]), 100 + i);
    }

    for(int64_t i = 0; i < 4; i++)
        delete callbacks[i];

}
