//
// Created by Matthew Nicoletti on 5/22/18.
//

#ifndef ANANTH_SAMPLING_SIMPLE_TABLEAUX_H
#define ANANTH_SAMPLING_SIMPLE_TABLEAUX_H

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "Tiling.h"


using namespace std;


class SimpleTableaux : public Tiling {
public:
    vector<double> vars;
    vector<int> lambda;
    vector<int> mu;
    int W;
    unsigned long h;
    double q;
    int m;
    int t;

    SimpleTableaux();

    SimpleTableaux(vector<int> &init, int n_init, int Width, unsigned long height, int tval, int mval);

    void set_to_top(int n_init, int tval, int mval);

    void set_to_bottom(int n_init, int tval, int mval);

    void set_params(vector<int>& init, int n_init, int tval, int mval);

    int get(int i, int j) override;

    void set(int i, int j, int val);

    double f(int i, int j);

    int c(int i, int j);

    int max_val(int i, int j);

    int min_val(int i, int j);

    double face_weight(int val, int i, int j);

    void step(double r, unsigned long seed) override;

    void transform(double r, int i, int j);

    void print() override;

    void print_f_vals() override;

    bool equals(Tiling& T) override;

    double max_diff(Tiling& T) override;

//    Tiling& operator=(const Tiling& other) override;
};


#endif //ANANTH_SAMPLING_SIMPLE_TABLEAUX_H
