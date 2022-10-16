#pragma once

#include "setup_tests.hh"
#include "delta_vector.hh"

TEST(DELTA_VECTOR, basic) {
    vector<int64_t> v = {0, 3, 6, 7, 12, 13, 14, 1000, 1000000};
    Delta_Vector dv(v);
    vector<int64_t> v2 = dv.get_values();
    ASSERT_EQ(v, v2);
}