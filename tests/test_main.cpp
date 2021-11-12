#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include <cassert>
#include <set>
#include <unordered_map>
#include <map>
#include <gtest/gtest.h>
#include "stdlib_printing.hh"
#include "globals.hh"
#include "setup_tests.hh"
#include <cassert>

#include "test_misc.hh"
#include "Kmer_tests.hh"
#include "test_coloring.hh"
#include "test_EM_sort.hh"
#include "test_pseudoalignment.hh"
#include "BOSS_tests.hh"
#include "test_extract_unitigs.hh"

int main(int argc, char **argv) {
    setup_tests(argc, argv);
    return RUN_ALL_TESTS();
}