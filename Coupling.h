//
// Created by Matthew Nicoletti on 6/22/18.
//
#include <iostream>
#include <vector>
#include "Tiling.h"

#ifndef ANANTH_SAMPLING_COUPLING_H
#define ANANTH_SAMPLING_COUPLING_H

using namespace std;


class Coupling {
public:
    Coupling(Tiling& bottom, Tiling& bot_copy, Tiling& top, Tiling& top_copy);
    Tiling& Top;
    Tiling& Bottom;
    Tiling& Bottom_copy;
    Tiling& Top_copy;
    Tiling& result;
    void cftp(int seed_of_seeds, int first_step_num);
    double check_dist();

};


#endif //ANANTH_SAMPLING_COUPLING_H
