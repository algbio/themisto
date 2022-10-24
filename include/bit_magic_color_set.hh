#pragma once

#include <iostream>
#include <vector>

#include <cstdint>
#include <cstring>

#include "bitmagic/bm64.h"
#include "bitmagic/bmserial.h"

class Bit_Magic_Color_Set {
    std::vector<unsigned char> vec;

public:
    Bit_Magic_Color_Set() {}

    Bit_Magic_Color_Set(const std::vector<unsigned char>& vec) : vec(vec) {}

    Bit_Magic_Color_Set(bm::bvector<> bv) {
        bv.optimize();

        bm::serializer<bm::bvector<> > bvs;
        bvs.byte_order_serialization(false);
        bvs.gap_length_serialization(false);

        bm::serializer<bm::bvector<>>::buffer sbuf;
        bvs.serialize(bv, sbuf);

        vec.resize(sbuf.size());
        std::memcpy(vec.data(), sbuf.data(), sbuf.size());
    }

    Bit_Magic_Color_Set(const std::vector<std::int64_t>& colors) {
        bm::bvector<> bv;

        for (const auto x : colors)
            bv.set(static_cast<std::uint64_t>(x));

        bv.optimize();

        bm::serializer<bm::bvector<> > bvs;
        bvs.byte_order_serialization(false);
        bvs.gap_length_serialization(false);

        bm::serializer<bm::bvector<>>::buffer sbuf;
        bvs.serialize(bv, sbuf);

        vec.resize(sbuf.size());
        std::memcpy(vec.data(), sbuf.data(), sbuf.size());
    }

    bool empty() const {
        return vec.empty();
    }

    std::int64_t size() const {
        if (empty())
            return 0;

        bm::bvector<> bv;
        bm::deserialize(bv, vec.data());

        return static_cast<std::int64_t>(bv.count());
    }

    std::int64_t size_in_bits() const {
        return static_cast<std::int64_t>(sizeof(std::vector<unsigned char>) + vec.size() * 8);
    }

    bool contains(const std::int64_t color) const {
        if (empty())
            return false;

        bm::bvector<> bv;
        bm::deserialize(bv, vec.data());

        return bv.get_bit(color);
    }

    void intersection(const Bit_Magic_Color_Set& other){
        if (empty()) return;
        if (other.empty()){
            *this = other;
            return;
        }

        bm::operation_deserializer<bm::bvector<>> od;
        bm::bvector<> bv;

        bm::deserialize(bv, vec.data());
        od.deserialize(bv, other.vec.data(), bm::set_AND);

        *this = Bit_Magic_Color_Set(bv);
    }

    void do_union(const Bit_Magic_Color_Set& other){
        if (empty()){
            *this = Bit_Magic_Color_Set(other.vec);
            return;
        }
        else if (other.empty()){
            *this = Bit_Magic_Color_Set(vec);
            return;
        } 

        bm::operation_deserializer<bm::bvector<>> od;
        bm::bvector<> bv;

        bm::deserialize(bv, vec.data());
        od.deserialize(bv, other.vec.data(), bm::set_OR);

        *this = Bit_Magic_Color_Set(bv);
    }

    std::vector<std::int64_t> get_colors_as_vector() const {
        if (empty())
            return std::vector<std::int64_t>();

        bm::bvector<> bv;
        bm::deserialize(bv, vec.data());

        return std::vector<std::int64_t>(bv.first(), bv.end());
    }

    std::int64_t serialize(std::ostream& os) const {
        std::size_t sz = vec.size();

        os.write(reinterpret_cast<char*>(&sz), sizeof(std::size_t));
        os.write(reinterpret_cast<const char*>(vec.data()), sz);

        return sizeof(std::size_t) + sz;
    }

    void load(std::istream& is) {
        std::size_t sz;
        is.read(reinterpret_cast<char*>(&sz), sizeof(std::size_t));

        vec.resize(sz);
        is.read(reinterpret_cast<char*>(vec.data()), sz);
    }
};
