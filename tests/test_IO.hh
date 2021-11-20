#pragma once

#include <gtest/gtest.h>
#include "../globals.hh"
#include "setup_tests.hh"
#include "test_tools.hh"
#include "input_reading.hh"
#include "ReadBatch.hh"

TEST(TEST_BUFFERED_IO, write_and_read){
    string filename = get_temp_file_manager().create_filename();

    LL n = 1e6 * 5; // Code below assumes that this is a multiple of ten
    vector<char> data(n);
    for(LL i = 0; i < n; i++) data[i] = rand() % 256;

    {
        Buffered_ofstream out(filename, ios::binary);
        
        out.set_buffer_capacity(123);
        
        for(LL i = 0; i < n; i += 10)
            out.write(data.data() + i, 10);
    } // End of scope should flush and close the stream

    vector<char> data2(n);
    
    Buffered_ifstream in(filename, ios::binary);
    in.set_buffer_capacity(114);

    for(LL i = 0; i < n; i += 10)
        in.read(data2.data() + i, 10);

    ASSERT_EQ(data, data2);
} 