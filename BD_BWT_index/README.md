Only supports inputs up to size 2^31 - 1 bytes

First compile sdsl-lite with:

```
git submodule init
git submodule update
cd sdsl-lite
sh ./install.sh
```

Then compile the rest with:

```
cmake -DCMAKE_BUILD_TYPE=Release .
make
```

or:
```
cmake -DCMAKE_BUILD_TYPE=Debug . 
make
```

Library files will be compiled into ./lib and
headers will be installed to ./include

To use the library, link against -lbdbwt -ldbwt -ldivsufsort64 -lsdsl
For example: 

```
g++ main.cpp -std=c++11 -L lib -I include -lbdbwt -ldbwt -ldivsufsort64 -lsdsl
```

The header include/BD_BWT_index.hh contains some documentation.

For example, below is code that iterates all right-maximal
nodes of the suffix link tree of the string mississippi

```
#include <iostream>
#include "include/BD_BWT_index.hh"
#include <string>
#include <set>

using namespace std;

set<char> alphabet;

void traverse(BD_BWT_index<>& index, Interval_pair I, string rev_label){
    for(char c : alphabet){
        Interval_pair I2 = index.left_extend(I,c);
        if(I2.forward.size() != 0 && index.is_right_maximal(I2)){
            rev_label += c;

            // Print the label of the current string
            string label(rev_label.rbegin(), rev_label.rend());
            cout << label << endl;

            // Recurse
            traverse(index, I2, rev_label);
            rev_label.pop_back(); // Undo extension
        }
    }
}

int main(){
    string text = "mississippi";
    for(auto c : text) alphabet.insert(c);
    BD_BWT_index<> index((uint8_t*)text.c_str());
    traverse(index, Interval_pair(0, text.size(), 0, text.size()), "");
}
```


