//
// Created by Matthew Nicoletti on 5/22/18.
//

#include "SimpleTableaux.h"
#include <random>
#include <cmath>

SimpleTableaux::SimpleTableaux() = default;

SimpleTableaux::SimpleTableaux(vector<int> &init, int n_init, int Width, unsigned long height, int tval, int mval) {
    m = mval;
    t = tval;
    n = n_init;
    W = Width;
    h = height;
    vector<int> tab((h + 2) * (W + 2));


    for (unsigned i = 0; i < height + 2; ++i) {
        for (unsigned j = 0; j < W + 2; ++j) {
            if ((j == 0) or (i == 0)) {
                tab[j + i * (W + 2)] = t * m * (n + j - i) - j;
            }
            else if (j == W + 1 or i == height + 1) {
                tab[j + i * (W + 2)] = -1;
            } else {
                tab[j + i * (W + 2)] = init[(j-1) + (i-1) * W];
            }
        }
    }
    tableaux = tab;
}


void SimpleTableaux::set_params(vector<int> &init, int n_init, int tval, int mval) {
    m = mval;
    t = tval;
    n = n_init;
    W = n;
    h = (unsigned long) n;
    vector<int> tab((h + 2) * (W + 2));


    for (unsigned i = 0; i < h + 2; ++i) {
        for (unsigned j = 0; j < W + 2; ++j) {
            if ((j == 0) or (i == 0)) {
                tab[j + i * (W + 2)] = t * m * (n + j - i) - j;
            }
            else if (j == W + 1 or i == h + 1) {
                tab[j + i * (W + 2)] = -1;
            } else {
                tab[j + i * (W + 2)] = init[(j-1) + (i-1) * W];
            }
        }
    }
    tableaux = tab;
}

void SimpleTableaux::set_to_top(int n_init, int tval, int mval) {
    m = mval;
    t = tval;
    n = n_init;

    tableaux = vector<int>((h + 2) * (W + 2));

    for (int i = 0; i < h + 2; ++i) {
        for (int j = 0; j < W + 2; ++j) {
            if ((j == 0) or (i == 0)) {
                tableaux[j + i * (W + 2)] = t * m * (n + j - i) - j;
            }
            else if (j == W + 1 or i == h + 1) {
                tableaux[j + i * (W + 2)] = -1;
            } else {
                tableaux[j + i * (W + 2)] = max_val(i, j);
            }
        }
    }
}

void SimpleTableaux::set_to_bottom(int n_init, int tval, int mval) {
    m = mval;
    t = tval;
    n = n_init;
//    W = n;
//    h = (unsigned long) n;
    tableaux = vector<int>((h + 2) * (W + 2));

    for (int i = h + 1; i >= 0; i--) {
        for (int j = W + 1; j >= 0; j--) {
            if ((j == 0) or (i == 0)) {
                tableaux[j + i * (W + 2)] = t * m * (n + j - i) - j;
            }
            else if (j == W + 1 or i == h + 1) {
                tableaux[j + i * (W + 2)] = -1;
            } else {
                tableaux[j + i * (W + 2)] = min_val(i, j);
            }
        }
    }
}

/* Tableaux's are 1-indexed (not 0-indexed) */
int SimpleTableaux::get(int i, int j) {
    return tableaux[i * (W + 2) + j];
}

void SimpleTableaux::set(int i, int j, int val) {
    tableaux[i * (W + 2) + j] = val;
}

/* If s = (i, j) is the position of a square which is not in the Tableaux, behavior is undefined.*/
double SimpleTableaux::f(int i, int j) {
    if (n + c(i, j) == 0) return t * m;
    return get(i, j) * 1.0 / (n + c(i, j));
}

int SimpleTableaux::c(int i, int j) {
    return j - i;
}

///* These correspond to Tableaux in LHT( >, > ) */
//int SimpleTableaux::max_val(int i, int j) {
//    int left, above;
//    double d = (n + c(i, j));
//    double e1 = .5 / (n + c(i - 1, j));
//    double e2 = .5 / (n + c(i, j - 1));
//    above = (int) floor(f(i - 1, j) * d - e1);
//    left = (int) floor(f(i, j - 1) * d - e2);
//    return max(min(left, above), 0);
//}
//
//int SimpleTableaux::min_val(int i, int j) {
//    int right, below;
//    double d = (n + c(i, j));
//    double e1 = .5 / (n + c(i + 1, j));
//    double e2 = .5 / (n + c(i, j + 1));
//    below = (int) ceil(f(i + 1, j) * d + e1);
//    right = (int) ceil(f(i, j + 1) * d + e2);
//    return max(max(below, right), 0);
//}

/* These correspond to Tableaux in LHT( >=, > ) */
int SimpleTableaux::max_val(int i, int j) {
    int left, above;
    double d = (n + c(i, j));
    double e = .5 / (n + c(i - 1, j));
    above = (int) floor(f(i - 1, j) * d - e);
    left = (int) floor(f(i, j - 1) * d);
    return min(left, above);
}

int SimpleTableaux::min_val(int i, int j) {
    int right, below;
    double d = (n + c(i, j));
    double e = .5 / (n + c(i + 1, j));
    below = (int) ceil(f(i + 1, j) * d + e);
    right = (int) ceil(f(i, j + 1) * d);
    return max(max(below, right), 0);
}

///* These correspond to Tableaux in LHT( >=, >= ) */
//int SimpleTableaux::max_val(int i, int j) {
//    int left, above;
//    double d = (n + c(i, j));
//    above = (int) floor(f(i - 1, j) * d);
//    left = (int) floor(f(i, j - 1) * d);
//    return min(left, above);
//}
//
//int SimpleTableaux::min_val(int i, int j) {
//    int right, below;
//    double d = (n + c(i, j));
//    below = (int) ceil(f(i + 1, j) * d);
//    right = (int) ceil(f(i, j + 1) * d);
//    return max(max(below, right), 0);
//}

double SimpleTableaux::face_weight(int val, int i, int j) {
    return 1;
}


void SimpleTableaux::step(double r, unsigned long seed) {
    int color = r < 0.5;
    int k = 0;
    srand(seed);
    for(unsigned i=1; i <= h; ++i) {
        for (unsigned j = 1; j <= W; ++j) {
            if ((i + j) % 2 == color) {
                transform(rand()*1.0/(RAND_MAX), i, j);
            }
        }
    }
}

void SimpleTableaux::transform(double r, int i, int j) {
    int min_v = min_val(i, j);
    int max_v = max_val(i, j);
    int new_val = (int) floor(r * (max_v - min_v + 1)) + min_v;
    set(i, j, new_val);
}

void SimpleTableaux::print() {
    cout << "{" << endl;
    for (unsigned i=0; i < h + 1; ++i) {
        cout << "{";
        for (unsigned j = 0; j < W + 1; ++j) {
            cout << get(i, j);
            if (j != W) {
                cout << ", ";
            }
        }
        cout << "}";
        if (i != h) cout << ",";
        cout << '\n';

    }
    cout << "}" << endl;
}

void SimpleTableaux::print_f_vals() {
    cout << "{" << endl;
    for (unsigned i = 0; i < h + 1; ++i) {
        cout << "{";
        for (unsigned j = 0; j < W + 1; ++j) {
            cout << f(i, j);
            if (j != W) {
                cout << ", ";
            }
        }
        cout << "}";
        if (i != h) cout << ",";
        cout << '\n';

    }
    cout << "}" << endl;
}

bool SimpleTableaux::equals(Tiling& T) {
    return tableaux == T.tableaux;
}

double SimpleTableaux::max_diff(Tiling& T) {
    double max_diff = 0.0;
    for (unsigned i=0; i < h + 1; ++i) {
        for (unsigned j = 0; j <= W; ++j) {
            max_diff = max( max_diff, abs(f(i, j) - T.get(i, j) * 1.0/(n + c(i, j))) );
        }
    }
    return max_diff;
}

//Tiling& SimpleTableaux::operator=(const Tiling& other) {
//    tableaux = other.tableaux;
//    return *this;
//}
