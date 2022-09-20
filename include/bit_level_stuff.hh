#pragma once

#include "globals.hh"

inline LL byte_to_int(char c){
    return static_cast<int>(*reinterpret_cast<unsigned char*>(&c));
}

inline LL parse_big_endian_LL(const char* A){
    LL x = 0;
    x |= byte_to_int(A[7]) << 0;
    x |= byte_to_int(A[6]) << 8;
    x |= byte_to_int(A[5]) << 16;
    x |= byte_to_int(A[4]) << 24;
    x |= byte_to_int(A[3]) << 32;
    x |= byte_to_int(A[2]) << 40;
    x |= byte_to_int(A[1]) << 48;
    x |= byte_to_int(A[0]) << 56;
    return x;
}

// Get the byte with index byte_idx in the big_endian representation of x
inline char get_byte(LL x, LL byte_idx){
    // How does castring int to char work?
    // Quote from the book The C++ Programming Language, by Bjarne Stroustrup:
    // > If the destination type is unsigned, the resulting value is simply as many bits from
    // > the source as will fit in the destination (high-order bits are thrown away if necessary).
    // > More precisely, the result is the least unsigned integer congruent to the source integer
    // > modulo 2 to the n-th, where n is the number of bits used to represent the unsigned type.
    // > If the destination type is signed, the value is unchanged if it can be represented in the
    // > destination type, otherwise, the value is implementation-defined.

    // Therefore, we can cast to an unsigned char and we should be ok. But the code that calls
    // this function needs a signed char, because ofstream::write takes only signed chars.
    // Therefore, to avoid implementation-defined behavior, we need to do a reinterpret-cast
    // to a signed value.
    unsigned char c = static_cast<unsigned char>((x >> (7 - byte_idx)*8) & 0xff);
    return *reinterpret_cast<char*>(&c);
}

inline void write_big_endian_LL(Buffered_ofstream<std::ofstream>& out, LL x){
    char c;

    c = get_byte(x,0); out.write(&c,1);
    c = get_byte(x,1); out.write(&c,1);
    c = get_byte(x,2); out.write(&c,1);
    c = get_byte(x,3); out.write(&c,1);
    c = get_byte(x,4); out.write(&c,1);
    c = get_byte(x,5); out.write(&c,1);
    c = get_byte(x,6); out.write(&c,1);
    c = get_byte(x,7); out.write(&c,1);
}


inline void write_big_endian_LL(char* buf, LL x){
    buf[0] = get_byte(x,0);
    buf[1] = get_byte(x,1);
    buf[2] = get_byte(x,2);
    buf[3] = get_byte(x,3);
    buf[4] = get_byte(x,4);
    buf[5] = get_byte(x,5);
    buf[6] = get_byte(x,6);
    buf[7] = get_byte(x,7);
}