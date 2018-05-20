//
// Created by Matthew Nicoletti on 5/16/18.
//
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#ifndef ANANTH_SAMPLING_TABLEAUX_H
#define ANANTH_SAMPLING_TABLEAUX_H

using namespace std;


class Tableaux {
    public:
    vector<double> vars;
    double u, v;
    vector<int> lambda;
    vector<int> mu;
    vector<int> tableaux;
    int length;
    int N;
    unsigned long height;
    int n;
    double q;

    Tableaux(vector<int>& l, vector<int>& m, vector<int>& init, int n_init, double u0, double v0, double qVal);

    int get(int i, int j);

    void set(int i, int j, int val);

    double f(int i, int j);

    int c(int i, int j);

    int max_val(int i, int j);

    int min_val(int i, int j);

    double face_weight(int val, int i, int j);

    int o(int m);

    double x(int i);

    double partition_fn(int i, int j);

};


#endif //ANANTH_SAMPLING_TABLEAUX_H
