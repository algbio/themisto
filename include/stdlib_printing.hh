#pragma once

#include <unordered_map>
#include <map>
#include <vector>
#include <set>
#include <utility>
#include <iostream>

using namespace std;

template <typename S, typename T>
ostream& operator<<(ostream& os, const unordered_map<S,T>& v){
    os << "[";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << it->first << ": " << it->second;
    }
    os << "]";
    return os;
}

template <typename S, typename T>
ostream& operator<<(ostream& os, const map<S,T>& v){
    os << "{";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << it->first << ": " << it->second;
    }
    os << "}";
    return os;
}

template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v){
    os << "[";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << *it;
    }
    os << "]";
    return os;
}

template <typename T>
ostream& operator<<(ostream& os, const set<T>& v){
    os << "[";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << *it;
    }
    os << "]";
    return os;
}

template <typename T>
ostream& operator<<(ostream& os, const multiset<T>& v){
    os << "[";
    for(auto it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) os << ", ";
        os << *it;
    }
    os << "]";
    return os;
}

template <typename S, typename T>
ostream& operator<<(ostream& os, const pair<S,T>& x){
    os << "(" << x.first << ", " << x.second << ")";
    return os;
}
