#include <iostream>
#include <vector>
#include <cmath>
#include "Tableaux.h"
#include <ctime>       /* time */

using namespace std;

void print(Tableaux& T);


int main() {
    unsigned int seed = 2;
    int n = 5;
    double u = .8;
    double v = .25;
    double q = .5;

    vector<int> lambda ({6, 6, 4, 3});
    vector<int> mu ({3, 1});
    vector<int> t ({
                           0, 0, 0, 9, 4, 3,
                           0, 5, 6, 4, 3, 1,
                           2, 2, 1, 0, 0, 0,
                           1, 0, 0, 0, 0, 0
                   });
    vector<double> x;

    Tableaux T(lambda, mu, t, n, u, v, q);
    print(T);
    T.step(seed);
//    random_device rd;
//    mt19937 mt(rd());
//    uniform_real_distribution<double> dist(0.0, 1.0);
//    double r = dist(mt);
//    T.transform(r, 3, 1);
    print(T);
    return 0;
}

void print(Tableaux& T) {
    for (unsigned i=0; i < T.height + 2; ++i) {
        for (unsigned j = 0; j < T.N + 2; ++j) {
            cout << T.get(i, j) << ' ';
        }
        cout << endl;
    }
}
