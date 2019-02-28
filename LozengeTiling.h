//
// Created by Matthew Nicoletti on 5/22/18.
//

#ifndef ANANTH_SAMPLING_LOZENGETILING_H
#define ANANTH_SAMPLING_LOZENGETILING_H

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "Tiling.h"

using namespace std;

class LozengeTiling : public Tiling {

public:
    LozengeTiling();

    LozengeTiling(int Nval, vector<int>& init, int n_init, int tval, int mval, int h);

    int N;
    int height;
    int m;
    int t;

    void set_params(int Nval, vector<int>& init, int n_init, int tval, int mval, int h);

    int get(int i, int j) override;

    void set(int i, int j, int val);

    int max_val(int i, int j);

    int min_val(int i, int j);

    void step(double r, unsigned long seed) override;

    void transform(double r, int i, int j);

    void print() override;

    bool equals(Tiling& T) override;

    double max_diff(Tiling& T) override;

    void print_f_vals() override;

//    Tiling& operator=(const Tiling& other) override;

};


#endif //ANANTH_SAMPLING_LOZENGETILING_H
