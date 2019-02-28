//
// Created by Matthew Nicoletti on 6/23/18.
//

#include <iostream>
#include <vector>


#ifndef ANANTH_SAMPLING_TILING_H
#define ANANTH_SAMPLING_TILING_H

using namespace std;

class Tiling {
public:
//    virtual Tiling& operator=(const Tiling& other) = 0;
    vector<int> tableaux;
    int n;
    virtual int get(int i, int j) = 0;
    virtual void step(double r, unsigned long seed) = 0;
    virtual void print() = 0;
    virtual bool equals(Tiling& T) = 0;
    virtual double max_diff(Tiling& T) = 0;
    virtual void print_f_vals() = 0;
};


#endif //ANANTH_SAMPLING_TILING_H
