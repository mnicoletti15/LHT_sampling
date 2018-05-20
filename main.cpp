#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "Tableaux.h"
#include <ctime>       /* time */

using namespace std;

void step(unsigned int seed, Tableaux& T);
void transform(double r2, Tableaux& T, int i, int j);
void print(Tableaux& T);


int main() {
    unsigned int seed = 1;
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
    step(seed, T);
    print(T);
    return 0;
}

void step(unsigned int seed, Tableaux& T) {
    srand(seed);

    double r = ((double) rand() / (RAND_MAX));
    int k = 0;
    for(unsigned i=1; i <= T.height; ++i) {
        for (unsigned j = 1; j <= T.N; ++j) {
            if ((i <= T.mu.size() and j <= T.mu[i - 1]) or (j > T.lambda[i - 1])) {
                continue;
            } else {
                k += 1;
                if (r <= (double) k / T.length) {
                    double r2 = ((double) rand() / (RAND_MAX));
                    transform(r2, T, i, j);
                    return;
                }
            }
        }
    }

}

void transform(double r2, Tableaux& T, int i, int j) {
    int k = T.min_val(i, j);
    double Z = T.partition_fn(i, j);
    cout << T.get(0, 0) << endl;
    cout << Z << endl;
    double c = T.face_weight(k, i, j)/Z;
    while (c < r2) {
        c += T.face_weight(++k, i, j)/Z;
    }
    T.set(i, j, k);
}

void print(Tableaux& T) {
    for (unsigned i=0; i < T.height + 2; ++i) {
        for (unsigned j = 0; j < T.N + 2; ++j) {
            cout << T.get(i, j) << ' ';
        }
        cout << endl;
    }
}
