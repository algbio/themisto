#pragma once

#include <iostream>
#include <gtest/gtest.h>
#include "globals.hh"
#include "sbwt/globals.hh"
#include "version.h"

class TestLogger{
    public:
    bool verbose = false;

    // This is to make std::endl compile and work with the logger
    TestLogger& operator<<(std::ostream& (*f)(std::ostream&)){
        if(verbose) f(std::cerr);
        return *this;
    }
};

template <typename T>
TestLogger& operator<<(TestLogger& L, const T& t){
    if(L.verbose) cerr << t;
    return L;
}

TestLogger logger; // Pipe things you want to print into this object with the '<<' operator

void enable_test_logging(){logger.verbose = true; }
void disable_test_logging(){logger.verbose = false; }

void setup_tests(int argc, char** argv){

    sbwt::write_log("Themisto-" + std::string(THEMISTO_BUILD_VERSION), sbwt::LogLevel::MAJOR);
    sbwt::write_log("Built at " + std::string(THEMISTO_BUILD_TIMESTAMP), sbwt::LogLevel::MAJOR);
    
    sbwt::write_log("Maximum k-mer length: " + std::to_string(MAX_KMER_LENGTH), sbwt::LogLevel::MAJOR);
    if(MAX_KMER_LENGTH != 255){
        throw std::runtime_error("Error: tests must be compiled with -DMAX_KMER_LENGTH=255 to cmake"); // 255 is kmer length, 255+1 is edgemer length
    }

    if(system("mkdir -p temp") != 0){
        cerr << "Error creating directory ./temp" << endl;
        exit(1);
    }

    if(system("mkdir -p test_data") != 0){
        cerr << "Error creating directory ./test_data" << endl;
        exit(1);
    }

    if(system("mkdir -p test_out") != 0){
        cerr << "Error creating directory ./test_out" << endl;
        exit(1);
    }
    
    bool verbose = false;
    for(int64_t i = 1; i < argc; i++)
        if(argv[i] == string("--verbose") || argv[i] == string("-v")) verbose = true;

    sbwt::get_temp_file_manager().set_dir("temp");

    verbose ? enable_test_logging() : disable_test_logging(); // test logger
    verbose ? sbwt::set_log_level(sbwt::LogLevel::DEBUG) : sbwt::set_log_level(sbwt::LogLevel::OFF); // main logger

    ::testing::InitGoogleTest(&argc, argv);

    srand(247829347);

}