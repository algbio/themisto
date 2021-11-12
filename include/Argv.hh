#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include <cassert>
#include <vector>

using namespace std;

class Argv{ // Class for turning a vector<string> into char**
private:

    // Forbid copying the class because it wont work right
    Argv(Argv const& other);
    Argv& operator=(Argv const& other);

public:

    char** array = NULL;
    int64_t size = 0;

    Argv(vector<string> v){
        array = (char**)malloc(sizeof(char*) * v.size());
        // Copy contents of v into array
        for(int64_t i = 0; i < v.size(); i++){
            char* s = (char*)malloc(sizeof(char) * (v[i].size() + 1)); // +1: space for '\0' at the end
            for(int64_t j = 0; j < v[i].size(); j++){
                s[j] = v[i][j]; // Can't use strcpy because s.c_str() is const
            }
            s[v[i].size()] = '\0';
            array[i] = s;
        }
        size = v.size();
    }

    ~Argv(){
        for(int64_t i = 0; i < size; i++) free(array[i]);
        free(array);
    }

};