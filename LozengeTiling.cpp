//
// Created by Matthew Nicoletti on 5/22/18.
//

#include "LozengeTiling.h"
#include <random>

using namespace std;

LozengeTiling::LozengeTiling() = default;

LozengeTiling::LozengeTiling(int Nval, vector<int> &init, int n_init, int tval, int mval, int h) {
    m = mval;
    t = tval;
    n = n_init;
    N = Nval;
    height = h;
    vector<int> tab((N + 2) * (height + 2));
    int W = N + 2;

    for (unsigned i = 0; i < height + 2; ++i) {
        for (unsigned j = 0; j < W; ++j) {
            if (j == 0 or i == 0) {
                tab[j + i * W] = t * m * n;
            }
            else if (j == W - 1 or i == height + 1) {
                tab[j + i * W] = -1;
            } else {
                tab[j + i * W] = init[(j-1) + (i-1) * N];
            }
        }
    }
    tableaux = tab;

}

void LozengeTiling::set_params(int Nval, vector<int> &init, int n_init, int tval, int mval, int h) {
    m = mval;
    t = tval;
    n = n_init;
    N = Nval;
    height = h;
    vector<int> tab((N + 2) * (height + 2));
    int W = N + 2;


    for (unsigned i = 0; i < height + 2; ++i) {
        for (unsigned j = 0; j < W; ++j) {
            if (j == 0 or i == 0) {
                tab[j + i * W] = t * m * n;
            }
            else if (j == W - 1 or i == height + 1) {
                tab[j + i * W] = -1;
            } else {
                tab[j + i * W] = init[(j-1) + (i-1) * N];
            }
        }
    }
    tableaux = tab;

}

/* Tableaux's are 1-indexed (not 0-indexed) */
int LozengeTiling::get(int i, int j) {
    return tableaux[i * (N + 2) + j];
}

void LozengeTiling::set(int i, int j, int val) {
    tableaux[i * (N + 2) + j] = val;
}

///* These correspond to Tableaux in ( >=, > ) */
//int LozengeTiling::max_val(int i, int j) {
//    int left, above;
//    left = get(i, j - 1);
//    above = get(i - 1, j);
//    return max(min(left, above - 1), 0);
//}
//
//int LozengeTiling::min_val(int i, int j) {
//    int right, below;
//    right = get(i, j + 1);
//    below = get(i + 1, j);
//    below = max(below, 0);
//    return max(max(below + 1, right), 0);
//}

/* this is with (>=, >=) */
int LozengeTiling::max_val(int i, int j) {
    int left, above;
    left = get(i, j - 1);
    above = get(i - 1, j);
    return min(left, above);
}

int LozengeTiling::min_val(int i, int j) {
    int right, below;
    right = get(i, j + 1);
    below = get(i + 1, j);
    below = max(below, 0);
    return max(right, below);
}

void LozengeTiling::step(double r, unsigned long seed) {
    srand(seed);
    int color = r < 0.5;
    for(unsigned i=1; i <= height; ++i) {
        for (unsigned j = 1; j <= N; ++j) {
            if ((i + j) % 2 == color) {
                transform(rand()*1.0/(RAND_MAX), i, j);
            }
        }
    }
}

void LozengeTiling::transform(double r, int i, int j) {
    int min_v = min_val(i, j);
    int max_v = max_val(i, j);
    int new_val = (int) floor(r * (max_v - min_v + 1)) + min_v;
    set(i, j, new_val);
}

void LozengeTiling::print() {
    cout << "{" << endl;
    for (unsigned i=0; i < height + 1; ++i) {
        cout << "{";
        for (unsigned j = 0; j <= N; ++j) {
            cout << get(i, j);
            if (j != N) {
                cout << ", ";
            }
        }
        cout << "}";
        if (i != height) cout << ",";
        cout << '\n';

    }
    cout << "}" << endl;
}

bool LozengeTiling::equals(Tiling& T) {
    return tableaux == T.tableaux;
}

double LozengeTiling::max_diff(Tiling& T) {
    int max_diff = 0;
    for (unsigned i=0; i < height + 1; ++i) {
        for (unsigned j = 0; j <= N; ++j) {
            max_diff = max(max_diff, abs(get(i, j) - T.get(i, j)));
        }
    }
    return max_diff * 1.0;
}

void LozengeTiling::print_f_vals() {
    cout << "no f vals for this tiling" << endl;
}

//
//Tiling& LozengeTiling::operator=(const Tiling& other) {
//    tableaux = other.tableaux;
//    return *this;
//}
