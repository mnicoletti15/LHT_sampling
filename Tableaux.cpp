//
// Created by Matthew Nicoletti on 5/16/18.
//

#include <limits>
#include "Tableaux.h"


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
                tab[j + i * W] = INFINITY;
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

double Tableaux::f(int i, int j) {
    return (double) get(i, j) / (n + c(i, j));
}

int Tableaux::c(int i, int j) {
    return j - i;
}

int Tableaux::max_val(int i, int j) {
    int left, above;
    double d = (n + c(i, j));
    double e = .5 / (n + c(i - 1, j));
    above = (int) floor(f(i - 1, j) * d - e);
    left = (int) floor(f(i, j - 1) * d);
    return min(left, above);
}

int Tableaux::min_val(int i, int j) {
    int right, below;
    double d = (n + c(i, j));
    double e = .5 / (n + c(i + 1, j));
    below = (int) ceil(f(i + 1, j) * d + e);
    right = (int) ceil(f(i, j + 1) * d);
    return max(below, right);
}

double Tableaux::face_weight(int val, int i, int j) {
    auto m = (int) floor(f(i, j));
    return x(val) * pow(u, m) * pow(v, o(m));
}

int Tableaux::o(int m) {
    return m % 2;
}

double Tableaux::x(int i) {
    return pow(q, i);
}

double geom_sum(double q, int lower, int upper) {
    return (pow(q, lower) - pow(q, upper + 1)) / (1 - q);
}

double Tableaux::partition_fn(int i, int j) {
    int mn = min_val(i, j);
    int mx = max_val(i, j);
    int M = n + c(i, j);
    auto a0 = (int) floor(mn * 1.0/M );
    auto a1 = (int) floor(mx * 1.0/M );
    int b0 = mn % M;
    int b1 = mx % M;
    auto up = (int) ceil((a1 - 1)/2.0);
    auto low = (int) ceil((a0 + 1)/2.0);

    double middle_chunk_even = geom_sum(pow(q, 2 * M)* pow(u, 2), low, up);
    double middle_chunk_odd = geom_sum(pow(q, 2 * M)* pow(u, 2), low - 1, up - 1);
    double middle_chunk = middle_chunk_even + v * pow(u, 2) * pow(q, M) * middle_chunk_odd;

    double result = 0.0;
    result += pow(u, a0) * pow(v, o(a0)) * pow(q, a0 * M) * geom_sum(q, b0, M-1);
    result += pow(u, a1) * pow(v, o(a1)) * pow(q, a1 * M) * geom_sum(q, 0, b1);
    result += middle_chunk * geom_sum(q, 0, M - 1);
    return result;
}
