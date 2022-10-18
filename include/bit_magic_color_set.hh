#pragma once

#include <iostream>
#include <vector>

#include <cstdint>

#include "bitmagic/bm64.h"
#include "bitmagic/bmserial.h"

class Bit_Magic_Color_Set {
    bm::bvector<> bv;
public:
    Bit_Magic_Color_Set() {}

    Bit_Magic_Color_Set(bm::bvector<> bv) : bv(bv) {}

    Bit_Magic_Color_Set(const std::vector<std::int64_t>& colors) {
        for (const auto x : colors)
            bv.set(static_cast<std::uint64_t>(x));

        bv.optimize();
    }

    bool empty() const {
        return bv.empty();
    }

    std::int64_t size() const {
        return static_cast<std::int64_t>(bv.count());
    }

    std::int64_t size_in_bits() const {
        bm::bvector<>::statistics st;
        bv.calc_stat(&st);

        return static_cast<std::int64_t>(st.memory_used * 8);
    }

    bool contains(const std::int64_t color) const {
        return bv.get_bit(color);
    }

    Bit_Magic_Color_Set intersection(const Bit_Magic_Color_Set& other) const {
        bm::bvector<> res;
        res.bit_and(bv, other.bv);

        return Bit_Magic_Color_Set(res);
    }

    Bit_Magic_Color_Set do_union(const Bit_Magic_Color_Set& other) const {
        bm::bvector<> res;
        res.bit_or(bv, other.bv);

        return Bit_Magic_Color_Set(res);
    }

    std::vector<std::int64_t> get_colors_as_vector() const {
        return std::vector<std::int64_t>(bv.first(), bv.end());
    }

    std::int64_t serialize(std::ostream& os) const {
        bm::serializer<bm::bvector<> > bvs;
        bvs.byte_order_serialization(false);
        bvs.gap_length_serialization(false);

        bm::serializer<bm::bvector<> >::buffer sbuf;
        bvs.serialize(bv, sbuf);
        std::size_t sz = sbuf.size();

        os.write(reinterpret_cast<char*>(&sz), sizeof(std::size_t));
        os.write(reinterpret_cast<char*>(sbuf.data()), sz);

        return sizeof(std::size_t) + sz;
    }

    void load(std::istream& is) {
        std::size_t sz;
        is.read(reinterpret_cast<char*>(&sz), sizeof(std::size_t));

        unsigned char* serialized_bytes = new unsigned char[sz];
        is.read(reinterpret_cast<char*>(serialized_bytes), sz);

        bm::deserialize(bv, serialized_bytes);
        bv.optimize();

        delete[] serialized_bytes;
    }
};
