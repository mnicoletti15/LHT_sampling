//
// Created by Matthew Nicoletti on 7/24/18.
//




//
// Created by Matthew Nicoletti on 5/16/18.
//
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "Tiling.h"
#ifndef ANANTH_SAMPLING_BOUNDEDTABLEAU_H
#define ANANTH_SAMPLING_BOUNDEDTABLEAU_H

using namespace std;


class BoundedTableaux : public Tiling {
public:
    double u, v;
    vector<int> lambda;
    vector<int> mu;
    int N;
    unsigned long h;
    double q;
    int m;
    int t;

    BoundedTableaux();

    BoundedTableaux(vector<int> &l, vector<int> &m, int n_init, double u0, double v0, double qVal);

    void set_params(vector<int>& l, vector<int>& m, vector<int>& init, int n_init, double u0, double v0, double qVal);

    void set_to_top(vector<int>& l, int n_init, int tval, int mval);

    void set_to_bottom(vector<int>& l, int n_init, int tval, int mval);

    int get(int i, int j) override;

    void set(int i, int j, int val);

    double f(int i, int j);

    int c(int i, int j);

    int max_val(int i, int j);

    int min_val(int i, int j);

    double face_weight(int val, int i, int j);

    int o(int m);

    double x(int i);

    void step(double r, unsigned long seed) override;

    void transform(double r, int i, int j);

    void print() override;

    bool equals(Tiling& T) override;

    double max_diff(Tiling& T) override;

    void print_f_vals() override;

};



#endif //ANANTH_SAMPLING_BOUNDEDTABLEAU_H
