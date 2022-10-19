#pragma once

#include "setup_tests.hh"
#include "delta_vector.hh"

template<typename delta_vector_t>
void test_delta_vector(){
    vector<int64_t> v = {0, 3, 6, 7, 12, 13, 14, 1000, 1000000};
    delta_vector_t dv(v);

    ASSERT_EQ(v, dv.get_values());
    ASSERT_FALSE(dv.empty());

    delta_vector_t dv2 = to_disk_and_back(dv);
    ASSERT_EQ(dv2.get_values(), dv.get_values());
    ASSERT_FALSE(dv2.empty());

    delta_vector_t dv3;
    ASSERT_TRUE(dv3.empty());

    delta_vector_t dv4 = {};
    ASSERT_TRUE(dv4.empty());

}

TEST(DELTA_VECTOR, basic) {
    test_delta_vector<Elias_Delta_Vector>();
    test_delta_vector<Fixed_Width_Delta_Vector>();
}