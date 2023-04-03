#pragma once

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include <cstdint>
#include <cstring>

#include <sdsl/bit_vectors.hpp>
#include <vector>
#include <cstdint>

#include "roaring/roaring.hh"
#include "roaring/roaring64map.hh"

class Roaring_Color_Set_View;

class Roaring_Color_Set {
    Roaring64Map roaring;

public:

    typedef Roaring_Color_Set_View view_t;

    Roaring_Color_Set() {}

    Roaring_Color_Set(const view_t& view); // This is defined at the .cpp file because the view is not yet defined here

    Roaring_Color_Set(const Roaring_Color_Set& r){
        *this = r;
    }

    Roaring_Color_Set(Roaring64Map r) : roaring(r) {}

    Roaring_Color_Set(const std::vector<std::int64_t>& colors) {
        for (const auto x : colors)
            roaring.add(static_cast<std::uint64_t>(x));

        roaring.runOptimize();
        roaring.shrinkToFit();
    }

    Roaring_Color_Set(const std::size_t n, const std::int32_t* colors) {
        roaring.addMany(n, reinterpret_cast<const std::uint32_t*>(colors));

        roaring.runOptimize();
        roaring.shrinkToFit();
    }

    Roaring_Color_Set(const std::size_t n, const std::int64_t* colors) {
        roaring.addMany(n, reinterpret_cast<const std::uint64_t*>(colors));

        roaring.runOptimize();
        roaring.shrinkToFit();
    }

    void add(const std::vector<std::int64_t>& colors) {
        for (const auto x : colors)
            roaring.add(static_cast<std::uint64_t>(x));

        roaring.runOptimize();
        roaring.shrinkToFit();
    }

    void add(const std::size_t n, const std::int32_t* colors) {
        roaring.addMany(n, reinterpret_cast<const std::uint32_t*>(colors));

        roaring.runOptimize();
        roaring.shrinkToFit();
    }

    void add(const std::size_t n, const std::int64_t* colors) {
        roaring.addMany(n, reinterpret_cast<const std::uint64_t*>(colors));

        roaring.runOptimize();
        roaring.shrinkToFit();
    }

    std::vector<std::int64_t> get_colors_as_vector() const {
        std::vector<std::int64_t> v(roaring.cardinality());
        roaring.toUint64Array(reinterpret_cast<std::uint64_t*>(v.data()));

        return v;
    }

    void push_colors_to_vector(std::vector<int64_t>& vec) const{
        // Todo: implement without get_colors_as_vector
        for(int64_t color : get_colors_as_vector())
            vec.push_back(color);
    }

    int64_t size() const {
        return roaring.cardinality();
    }

    bool empty() const{
        return size() == 0;
    }


    int64_t size_in_bits() const {
        return roaring.getSizeInBytes(false) * 8;
    }

    bool contains(const std::int64_t n) const {
        return roaring.contains((uint64_t)n);
    }

    void intersection(const Roaring_Color_Set& c) {
        *this = Roaring_Color_Set(roaring & c.roaring);
    }

    // union is a reserved word in C++ so this function is called do_union
    void do_union(const Roaring_Color_Set& c) {
        *this = Roaring_Color_Set(roaring | c.roaring);
    }

    int64_t serialize(std::ostream& os) const {
        std::size_t expected_size = roaring.getSizeInBytes(false);
        char* serialized_bytes = new char[expected_size];

        roaring.write(serialized_bytes, false);
        os.write(reinterpret_cast<char*>(&expected_size), sizeof(std::size_t));
        os.write(serialized_bytes, expected_size);
        delete[] serialized_bytes;

        return sizeof(std::size_t) + expected_size;
    }

    void load(std::istream& is) {
        std::size_t n;
        is.read(reinterpret_cast<char*>(&n), sizeof(std::size_t));

        char* serialized_bytes = new char[n];
        is.read(serialized_bytes, n);
        roaring = Roaring64Map::read(serialized_bytes, false);
        delete[] serialized_bytes;
    }
};

class Roaring_Color_Set_View {

    public:

    const Roaring_Color_Set* ptr; // Non-owning pointer

    Roaring_Color_Set_View(const Roaring_Color_Set& c) : ptr(&c) {}

    bool empty() const{
        return ptr->empty();
    }

    int64_t size() const{
        return ptr->size();
    }

    int64_t size_in_bits() const{
        return ptr->size_in_bits();
    }

    bool contains(int64_t color) const{
        return ptr->contains(color);
    }

    std::vector<int64_t> get_colors_as_vector() const{
        return ptr->get_colors_as_vector();
    }

    void push_colors_to_vector(std::vector<int64_t>& vec) const{
        return ptr->push_colors_to_vector(vec);
    }

};