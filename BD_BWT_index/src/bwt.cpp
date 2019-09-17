/*#include "bwt.hh"
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include "divsufsort64.h"

extern "C" { 
#include "dbwt.h"
}


using namespace std;

uint8_t* bwt_dbwt(uint8_t* text, int64_t length, uint8_t end_char){

    int64_t n = length;
    uint8_t* bwt = (uint8_t*) malloc((sizeof(uint8_t) * (n+2)));
    for(int64_t i = 0; i < length; i++)
        bwt[i] = text[i];
    bwt[n] = 0;
    bwt[n+1] = 0;
    
    int64_t last = divbwt64(bwt, bwt, NULL, n);
    for(int64_t i = n; i > last; i--) 
        bwt[i] = bwt[i-1];
    bwt[last] = end_char;
    
    return bwt;
}
*/