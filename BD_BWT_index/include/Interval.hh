#ifndef INTERVAL_HH
#define INTERVAL_HH

#include <string>
#include <sstream>
#include <utility>

class Interval{
public:
    int64_t left, right;
    Interval(int64_t left, int64_t right) : left(left), right(right) {}
    Interval() {}
    friend bool operator== (const Interval& x, const Interval& y);
    friend bool operator!= (const Interval& x, const Interval& y);
    friend bool operator<(const Interval& x, const Interval& y);
    int64_t size() const{
        return right - left + 1;
    }
    std::string toString(){
        std::stringstream ss;
        ss << "[" << left << ", " << right << "]";
        return ss.str();
    }
};

bool operator== (const Interval& x, const Interval& y){
    return x.left == y.left && x.right == y.right;
}

bool operator!= (const Interval& x, const Interval& y){
    return !(x == y);
}

bool operator<(const Interval& x, const Interval& y){
    return std::make_pair(x.left,x.right) < std::make_pair(y.left, y.right);
}

class Interval_pair{
public:
    Interval forward, reverse;
    Interval_pair(Interval forward, Interval reverse) : forward(forward), reverse(reverse) {}
    Interval_pair(int64_t left1, int64_t right1, int64_t left2, int64_t right2) :
        forward(left1,right1), reverse(left2, right2) {}
    Interval_pair() {}
    std::string toString(){
        return "(" + forward.toString() + ", " + reverse.toString() + ")";
    }
    friend bool operator== (const Interval_pair& x, const Interval_pair& y);
    friend bool operator!= (const Interval& x, const Interval& y);
    friend bool operator<(const Interval& x, const Interval& y);
};

bool operator== (const Interval_pair& x, const Interval_pair& y){
    return x.forward == y.forward && x.reverse == y.reverse;
}

bool operator!= (const Interval_pair& x, const Interval_pair& y){
    return !(x == y);
}

bool operator< (const Interval_pair& x, const Interval_pair& y){
    return std::make_pair(x.forward, x.reverse) < std::make_pair(y.forward, y.reverse);
}

#endif