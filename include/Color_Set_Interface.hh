#pragma once
#include <vector>

template <typename T>
concept Color_Set_Interface = requires(T& t){
    { t.empty() } -> std::same_as<bool>; // This should take constant time
    { t.size() } -> std::same_as<int64_t>; // This may take linear time
    { t.size_in_bits() } -> std::same_as<int64_t>;
    { t.contains(int64_t()) } -> std::same_as<bool>; // This may take linear time
    { t.intersection(t) };
    { t.do_union(t) };
    //{ t.serialize(os) } -> std::same_as<int64_t>; // Serialization not required for individual color sets any more
    //{ t.load(is) } -> std::same_as<void>;  // Loading not required for individual color sets any more
    { t.get_colors_as_vector() } -> std::same_as<std::vector<int64_t>>;
    requires std::constructible_from<T, std::vector<int64_t>>;
    requires std::default_initializable<T>;
};
