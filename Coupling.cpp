//
// Created by Matthew Nicoletti on 6/22/18.
//

#include "Coupling.h"
#include "Tiling.h"
#include <random>

using namespace std;

Coupling::Coupling(Tiling& bottom, Tiling& bot_copy, Tiling& top,  Tiling& top_copy) :
    Bottom(bottom),
    Bottom_copy(bot_copy),
    Top(top),
    Top_copy(top_copy),
    result(bot_copy)
{

}

void Coupling::cftp(int seed_of_seeds, int first_step_num) {
    mt19937 mt_rand(seed_of_seeds);
    vector<long> steps; vector<long> seeds;
    steps.push_back(first_step_num); seeds.push_back(mt_rand());

    bool coupled = false;
    while (true) {

        for (int i = 0; i < steps.size(); i++) {
            auto real_rand = bind(std::uniform_real_distribution<double>(0,1),
                                       mt19937(seeds[i]));

            for (int j = 0; j < steps[i]; j++) {
                double r= real_rand();
                unsigned long seed = rand();
                Top.step(r, seed);
                Bottom.step(r, seed);
            }
        }
//        coupled = Top.equals(Bottom);

        cout << Top.max_diff(Bottom) << endl;
        coupled = Top.max_diff(Bottom) < 2;

        if (coupled) {
            result = Top;
            break;
        }

        steps.insert(steps.begin(),steps[steps.size()-1]*2);
        seeds.insert(seeds.begin(),mt_rand());
        Top = Top_copy;
        Bottom = Bottom_copy;
    }
}

double Coupling::check_dist() {
    int s = 0;
    for (int i = 0; i < 10000; i++) {

    }
    return s * 1.0 / 10000.0;
}