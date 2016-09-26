//
// Created by Hermes235 on 26.9.2016.
//
#include <iostream>
#include <cstdlib>

using namespace std;

int main(){

    cout << "topology: 2 4 1" << endl;
    for (int i = 2000; i >= 0; --i) {
        int n1 = (int)(2.0 * rand() / double(RAND_MAX));
        int n2 = (int)(2.0 * rand() / double(RAND_MAX));
        int t = !(n1 && n2);

        cout << "in: " << n1 << ".0 " << n2 << ".0 " << endl;
        cout << "out: " << t << ".0" << endl;
    }
}

