#pragma once

#include "setup_tests.hh"
#include <gtest/gtest.h>
#include "globals.hh"
#include "sbwt/globals.hh"
#include "Sparse_Uint_Array.hh"

using namespace sbwt;

TEST(TEST_SPARSE_UINT_ARRAY, random_test){
    LL length = 1000;
    LL max_value = 9; // power of 2 plus 1 in case that is a special case
    Sparse_Uint_Array_Builder builder(length, 2048, 3);
    
    LL NOTFOUND = 1e9;
    vector<LL> reference(length, NOTFOUND);
    srand(12345);

    // Set 500 random values
    for(LL i = 0; i < 500; i++){
        
        LL index = rand() % length;
        LL value = rand() % (max_value+1);

        builder.add(index,value);

        reference[index] = min(reference[index], value);
        // ^ The comment in builder.add says that if multiple values go to the same index, the smallest values is kept
    }

    Sparse_Uint_Array A = builder.finish();

    // Check
    for(LL i = 0; i < length; i++){
        logger << i << " " << reference[i] << endl;
        if(reference[i] == NOTFOUND)
            ASSERT_EQ(A.get(i), -1);
        else
            ASSERT_EQ(A.get(i), reference[i]);
    }
}