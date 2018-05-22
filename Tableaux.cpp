//
// Created by Matthew Nicoletti on 5/16/18.
//

#include "Tableaux.h"
#include <random>

Tableaux::Tableaux(vector<int> &l, vector<int> &m, vector<int> &init, int n_init, double u0, double v0, double qVal) {
    q = qVal;
    u = u0;
    v = v0;
    n = n_init;
    lambda = l;
    mu = m;
    N = lambda.at(0);
    height = lambda.size();
    vector<int> tab((N + 2) * (height + 2));
    pair<vector<int>, int> result;
    int L = 0;
    int W = N + 2;


    for (unsigned i = 0; i < height + 2; ++i) {
        for (unsigned j = 0; j < W; ++j) {
            if ((j == 0) or (i == 0) or (i <= mu.size() and j <= mu[i - 1])) {
                tab[j + i * W] = -2;
            }
            else if ((j == W) or (i == height + 1) or (j > lambda[i - 1])) {
                tab[j + i * W] = -1;
            } else {
                tab[j + i * W] = init[(j-1) + (i-1) * N];
                L += 1;
            }
        }
    }
    tableaux = tab;
    length = L;

}

/* Tableaux's are 1-indexed (not 0-indexed) */
int Tableaux::get(int i, int j) {
    return tableaux[i * (N + 2) + j];
}

void Tableaux::set(int i, int j, int val) {
    tableaux[i * (N + 2) + j] = val;
}

/* If s = (i, j) is the position of a square which is not in the Tableaux, behavior is undefined.*/
double Tableaux::f(int i, int j) {
    return (double) get(i, j) * 1.0 / (n + c(i, j));
}

int Tableaux::c(int i, int j) {
    return j - i;
}

/* -2 is the integer which encodes +Infinity */
int Tableaux::max_val(int i, int j) {
    int left, above;
    double d = (n + c(i, j));
    double e = .5 / (n + c(i - 1, j));
    above = (int) floor(f(i - 1, j) * d - e);
    left = (int) floor(f(i, j - 1) * d);
    if (get(i - 1, j) == -2 and get(i, j - 1) == -2) return -2;
    else if (get(i - 1, j) == -2) return left;
    else if (get(i, j - 1) == -2) return above;
    else return min(left, above);
}

int Tableaux::min_val(int i, int j) {
    int right, below;
    double d = (n + c(i, j));
    double e = .5 / (n + c(i + 1, j));
    below = (int) ceil(f(i + 1, j) * d + e);
    right = (int) ceil(f(i, j + 1) * d);
    return max(max(below, right), 0);
}

double Tableaux::face_weight(int val, int i, int j) {
    int m = val / (n + c(i, j));
    return x(val) * pow(u, m) * pow(v, o(m));
}

int Tableaux::o(int m) {
    return m % 2;
}

double Tableaux::x(int i) {
    return pow(q, i);
}

double geom_sum(double q, int lower, int upper) {
    if (upper <= -2) {
        return pow(q, lower)/(1 - q);
    }
    if (upper < lower) {
        return 0;
    }
    else {
        return (pow(q, lower) - pow(q, upper + 1)) / (1 - q);
    }
}

double Tableaux::partition_fn(int i, int j) {
    int mn = min_val(i, j);
    int mx = max_val(i, j);

    int M = n + c(i, j);
    auto a0 = mn / M;
    auto a1 = mx / M;
    int b0 = mn % M;
    int b1 = mx % M;
    double result = 0.0;

//    if (mx == -2) {
//        result += pow(u, a0) * pow(v, o(a0)) * pow(q, a0 * M) * geom_sum(q, b0, M-1);
//    }


    auto up = (int) ceil((a1 - 1)/2.0);
    auto low = (int) ceil((a0 + 1)/2.0);

    if (mx == -2) {
        up = -2;
    }
//    cout << "a0: " << a0 << " a1: " << a1 << endl;
//    cout << "b0: " << b0 << " b1: " << b1 << endl;
    double middle_chunk_even = geom_sum(pow(q, 2 * M)* pow(u, 2), low, up);
    double middle_chunk_odd = geom_sum(pow(q, 2 * M)* pow(u, 2), low - 1, up - 1);
    double middle_chunk = middle_chunk_even + v * u * pow(q, M) * middle_chunk_odd;

    int first_upper_lim = M - 1;
    if (mx != -2 and a0 == a1) {
        first_upper_lim = b1;
    }
    result += pow(u, a0) * pow(v, o(a0)) * pow(q, a0 * M) * geom_sum(q, b0, first_upper_lim);
    if (mx != -2 and a1 > a0) {
        result += pow(u, a1) * pow(v, o(a1)) * pow(q, a1 * M) * geom_sum(q, 0, b1);
    }
    result += middle_chunk * geom_sum(q, 0, M - 1);
    return result;
}

void Tableaux::step(unsigned int seed) {
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);
    double r = dist(mt);
    int color = r < 0.5;

    int k = 0;
    for(unsigned i=1; i <= height; ++i) {
        for (unsigned j = 1; j <= N; ++j) {
            if ((i <= mu.size() and j <= mu[i - 1]) or (j > lambda[i - 1])) {
                continue;
            } else {
                k += 1;
                if ((i + j) % 2 == color) {
                    double r2 = dist(mt);
                    transform(r2, i, j);
                }
            }
        }
    }

}

void Tableaux::transform(double r, int i, int j) {
    int k = min_val(i, j);
    double Z = partition_fn(i, j);
    double c = face_weight(k, i, j)/Z;
    double s = 0;
    cout << "setting face " << i << ' ' << j << endl;
    while (c < r) {
        c += face_weight(++k, i, j)/Z;
    }
    cout << "done, set face " << i << ' ' << j << endl;
    set(i, j, k);
}
