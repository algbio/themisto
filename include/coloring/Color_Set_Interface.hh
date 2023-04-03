#pragma once
#include <vector>

template <typename T>
concept Color_Set_Interface = requires(T& t, typename T::view_t& view, std::vector<int64_t>& vec){
    { t.empty() } -> std::same_as<bool>; // This should take constant time
    { t.size() } -> std::same_as<int64_t>; // This may take linear time
    { t.size_in_bits() } -> std::same_as<int64_t>;
    { t.contains(int64_t()) } -> std::same_as<bool>; // This may take linear time
    { t.intersection(t) };
    { t.do_union(t) };
    { t.get_colors_as_vector() } -> std::same_as<std::vector<int64_t>>;
    { t.push_colors_to_vector(vec)};
    requires std::constructible_from<T, std::vector<int64_t>>;
    requires std::default_initializable<T>;

    // Below are the requirements for the associated view class T::view_t
    { view.empty() } -> std::same_as<bool>; // This should take constant time
    { view.size() } -> std::same_as<int64_t>; // This may take linear time
    { view.size_in_bits() } -> std::same_as<int64_t>;
    { view.contains(int64_t()) } -> std::same_as<bool>; // This may take linear time
    { view.get_colors_as_vector() } -> std::same_as<std::vector<int64_t>>;
    { view.push_colors_to_vector(vec)};
    requires std::constructible_from<typename T::view_t, T>;
    requires std::constructible_from<T, typename T::view_t>;
    // No intersection and union required for a view because a view is supposed to be immutable.
};
