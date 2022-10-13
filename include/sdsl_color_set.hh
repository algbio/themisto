#include "sdsl/bit_vectors.hpp"

class SDSL_Color_Set {
    sdsl::hyb_vector<> data;

public:

    SDSL_Color_Set(){}

    SDSL_Color_Set(const vector<std::int64_t>& colors) {
        int64_t max_color = *std::max_element(colors.begin(), colors.end());
        sdsl::bit_vector vec(max_color+1, 0);
        for(int64_t x : colors) vec[x] = 1;
        data = sdsl::hyb_vector<>(vec);
    }


    std::vector<std::uint64_t> get_colors_as_vector() const {
        std::vector<std::uint64_t> vec;
        for(int64_t i = 0; i < data.size(); i++){
            if(data[i]) vec.push_back(i);
        }
    }

    std::size_t size() const {
        return get_colors_as_vector().size(); // TODO SLOW
    }


    std::size_t size_in_bits() const {
        return sdsl::size_in_bytes(data) * 8;
    }

    bool contains(const std::uint64_t n) const {
        return data[n];
    }

    SDSL_Color_Set intersection(const SDSL_Color_Set& c) const {
        // TODO: SLOW
        vector<std::int64_t> result;
        int64_t n = min(data.size(), c.data.size());
        for(int64_t i = 0; i < n; i++){
            if(data[i] && c.data[i])
                result.push_back(i);
        }
        return SDSL_Color_Set(result);
    }

    // union is a reserved word in C++ so this function is called do_union
    SDSL_Color_Set do_union(const SDSL_Color_Set& c) const {
        // TODO: SUPER SLOW
        set<std::int64_t> result;
        for(std::int64_t x : get_colors_as_vector()){
            result.insert(x);
        }
        for(std::int64_t x : c.get_colors_as_vector()){
            result.insert(x);
        }

        vector<std::int64_t> vec(result.begin(), result.end());
        return SDSL_Color_Set(vec);
    }

    std::size_t serialize(std::ostream& os) const {
        return data.serialize(os);
    }

    std::size_t load(std::ifstream& is) {
        data.load(is);
        return sdsl::size_in_bytes(data);
    }
};