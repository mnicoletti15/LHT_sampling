//
// Created by Matthew Nicoletti on 5/16/18.
//

#include "BoundedTableaux.h"
#include <random>

BoundedTableaux::BoundedTableaux() = default;

BoundedTableaux::BoundedTableaux(vector<int> &l, vector<int> &m, int n_init, double u0, double v0, double qVal) {
    q = qVal;
    u = u0;
    v = v0;
    n = n_init;
    lambda = l;
    mu = m;
    N = lambda[1];
    h = lambda.size()-1;
    vector<int> tab((N + 2) * (h + 2));
    int W = N + 2;

    for (unsigned i = 0; i < h + 2; ++i) {
        for (unsigned j = 0; j < W; ++j) {
            if ((j == 0) or (i == 0) or (i <= mu.size() and j <= mu[i - 1])) {
                tab[j + i * W] = n * (n + j - i) - j;
            }
            else if ((j == W) or (i == h + 1) or (j > lambda[i])) {
                tab[j + i * W] = -1;
            } else {
                tab[j + i * W] = n * (n + j - i) - j;
            }
        }
    }
    tableaux = tab;

}

void BoundedTableaux::set_params(vector<int> &l, vector<int> &m, vector<int>& init, int n_init, double u0, double v0, double qVal) {
    q = qVal;
    u = u0;
    v = v0;
    n = n_init;
    lambda = l;
    mu = m;
    N = lambda[1];
    h = lambda.size()-1;
    vector<int> tab((N + 2) * (h + 2));
    int L = 0;
    int W = N + 2;


    for (unsigned i = 0; i < h + 2; ++i) {
        for (unsigned j = 0; j < W; ++j) {
            if ((j == 0) or (i == 0) or (i <= mu.size() and j <= mu[i - 1])) {
                tab[j + i * W] = n * (n + j - i) - j;
            }
            else if ((j == W) or (i == h + 1) or (j > lambda[i])) {
                tab[j + i * W] = -1;
            } else {
                tab[j + i * W] = init[(i-1) + N * (j-1)];
            }
        }
    }
    tableaux = tab;

}

void BoundedTableaux::set_to_top(vector<int> &l, int n_init, int tval, int mval) {
    n = n_init;
    lambda = l;
    N = lambda[1];
    h = lambda.size()-1;
    t = tval;
    m = mval;
    tableaux = vector<int>((h + 2) * (N + 2));
    int W = N;


    for (unsigned i = 0; i < h + 2; ++i) {
        for (unsigned j = 0; j < W + 2; ++j) {
            if ((j == 0) or (i == 0) or (i <= mu.size() and j <= mu[i - 1])) {
                tableaux[j + i * (W+2)] = t * m * (n + j - i) - j;
            }
            else if ((j == W) or (i == h + 1) or (j > lambda[i])) {
                tableaux[j + i * (W+2)] = -1;
            } else {
                tableaux[j + i * (W+2)] = t * m * (n + j - i) - j;
            }
        }
    }

}

void BoundedTableaux::set_to_bottom(vector<int> &l, int n_init, int tval, int mval) {
    n = n_init;
    lambda = l;
    N = lambda[1];
    h = lambda.size()-1;
    t = tval;
    m = mval;
    tableaux = vector<int>((h + 2) * (N + 2));
    int W = N;


    for (int i = h + 1; i >= 0; i--) {
        for (int j = W + 1; j >= 0; j--) {
            if ((j == 0) or (i == 0)) {
                tableaux[j + i * (W + 2)] = t * m * (n + j - i) - j;
            }
            else if (i == h + 1 or j > lambda[i]) {
                tableaux[j + i * (W + 2)] = -1;
            } else {
                tableaux[j + i * (W + 2)] = min_val(i, j);
            }
        }
    }

}

/* Tableaux's are 1-indexed (not 0-indexed) */
int BoundedTableaux::get(int i, int j) {
    return tableaux[i * (N + 2) + j];
}

void BoundedTableaux::set(int i, int j, int val) {
    tableaux[i * (N + 2) + j] = val;
}

/* If s = (i, j) is the position of a square which is not in the Tableaux, behavior is undefined.*/
double BoundedTableaux::f(int i, int j) {
    if (n + c(i, j) == 0) return t*m;
    return (double) get(i, j) * 1.0 / (n + c(i, j));
}

int BoundedTableaux::c(int i, int j) {
    return j - i;
}

int BoundedTableaux::max_val(int i, int j) {
    int left, above;
    double d = (n + c(i, j));
    double e = .5 / (n + c(i - 1, j));
    above = (int) floor(f(i - 1, j) * d - e);
    left = (int) floor(f(i, j - 1) * d);
    return min(left, above);
}

int BoundedTableaux::min_val(int i, int j) {
    int right, below;
    double d = (n + c(i, j));
    double e = .5 / (n + c(i + 1, j));
    below = (int) ceil(f(i + 1, j) * d + e);
    right = (int) ceil(f(i, j + 1) * d);
    return max(max(below, right), 0);
}


double BoundedTableaux::face_weight(int val, int i, int j) {
    int m = val / (n + c(i, j));
    return x(val) * pow(u, m) * pow(v, o(m));
}

int BoundedTableaux::o(int m) {
    return m % 2;
}

double BoundedTableaux::x(int i) {
    return pow(q, i);
}


void BoundedTableaux::step(double r, unsigned long seed) {
    int color = r < 0.5;
    int k = 0;
    srand(seed);
    for(unsigned i=1; i <= h; ++i) {
        for (unsigned j = 1; j <= lambda[i]; ++j) {
            if (i <= mu.size() and j <= mu[i - 1]) {
                continue;
            } else {
                k += 1;
                if ((i + j) % 2 == color) {
                    transform(rand() * 1.0/RAND_MAX, i, j);
                }
            }
        }
    }
}

void BoundedTableaux::transform(double r, int i, int j) {
    int min_v = min_val(i, j);
    int max_v = max_val(i, j);
    int new_val = (int) floor(r * (max_v - min_v + 1)) + min_v;
    set(i, j, new_val);
}

void BoundedTableaux::print() {
    cout << "{" << '\n';
    for (unsigned i=0; i < h + 2; ++i) {
        cout << "{";
        for (unsigned j = 0; j < N + 2; ++j) {
            cout << get(i, j) << ", ";
        }
        cout << "}" << endl;
    }
    cout << "}" << endl;
}

bool BoundedTableaux::equals(Tiling& T) {
    return tableaux == T.tableaux;
}

double BoundedTableaux::max_diff(Tiling& T) {
    double max_diff = 0.0;
    for (unsigned i=1; i < h + 1; ++i) {
        for (unsigned j = 1; j <= lambda[i]; ++j) {
            max_diff = max(max_diff, abs(get(i, j)*1.0/(n + j - i) - T.get(i, j)*1.0/(n + j - i)));
        }
    }
    return max_diff * 1.0;
}

void BoundedTableaux::print_f_vals() {
    cout << "{" << '\n';
    for (unsigned i=0; i < h + 1; ++i) {
        cout << "{";
        for (unsigned j = 0; j <= lambda[i]; ++j) {
            cout << f(i, j);
            if (j != lambda[i]) {
                cout << ", ";
            }
        }
        cout << "}";
        if (i != h) {
            cout << ", " << endl;
        }
    }
    cout << "}" << endl;
}

//Tiling& Tableaux::operator=(const Tiling& other) {
//    tableaux = other.tableaux;
//    return *this;
//}
