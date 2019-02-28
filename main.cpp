#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include "Tableaux.h"
#include "LozengeTiling.h"
#include "SimpleTableaux.h"
#include "Coupling.h"
#include "BoundedTableaux.h"
#include <ctime>       /* time */

using namespace std;

void print(Tableaux& T);
void print_tiling(LozengeTiling& T);
void run_LHT();
void print_tiling_derivative(LozengeTiling& T);
void run_Lozenge();
void run_simple();
void run_simple_coupling();
void run_Lozenge_coupling();
void run_bounded();
void run_one_path();
void run_skew_shape();
void run_skew_shape2();
void run_nonintersecting_coupling();
void run_small_t();
void cusps();
void rectangle_small_t();

int main() {

//    run_simple();
//    run_simple_coupling();
//    run_Lozenge();
//    run_bounded();
//    run_one_path();
//    run_skew_shape();
//    run_Lozenge_coupling();
//    run_nonintersecting_coupling();
    run_small_t();
//    run_skew_shape2();
//    rectangle_small_t();
//    cusps();

    return 0;
}

void rectangle_small_t() {
//    unsigned int seed = 21034;
    unsigned int seed = 2103044;
    int N = 100;
    int m = 1;
    int n = N * m;
    int W = N * m;
    unsigned long h = (unsigned long) N * m;
    srand(seed);

//    SimpleTableaux Top, Bottom, Top_copy, Bottom_copy;
//    Top.W = W; Top.h = h;
//    Bottom.W = W; Bottom.h = h;
//    Top.set_to_top(n, 100, m);
//    Bottom.set_to_bottom(n, 100, m);
//    Top_copy = Top;
//    Bottom_copy = Bottom;

    auto real_rand = bind(std::uniform_real_distribution<double>(0,1),
                          mt19937(seed));

    SimpleTableaux Top;
    Top.W = W; Top.h = h;
//    for (int t = 2; t <= 4; t++) {
//        Top.set_to_top(n, t, m);
//        for (int i = 0; i < pow(2, 16); i ++) {
//            Top.step(real_rand(), rand());
//
//        }
//        Top.print_f_vals();
//    }
    Top.set_to_top(n, 4, m);
    for (int i = 0; i < pow(2, 20); i ++) {
        Top.step(real_rand(), rand());
    }
    Top.print_f_vals();

}

void run_small_t() {
    unsigned int seed = 329392;
    int N = 60;
//    int t = 1;
    int t = 4;
    int m = 1;
    int n = N * m;
    int W = N * m;
    unsigned long h = (unsigned long) N * m;
    srand(seed);
    vector<int> lambda ((unsigned long) 2 * n + 1);
    lambda[0] = 8 * n;
    for (int i = 1; i <= n; i++) lambda[i] = 8*n - i;
    for (int i = n+1; i <= 2*n; i++) lambda[i] = 2*n - i + 1;
    vector<int> mu (0);
//    unsigned int seed = 21034;
//    int N = 20;
//    int t = 2;
//    int m = 2;
//    int n = N * m;
//    int W = N * m;
//    unsigned long h = (unsigned long) N * m;
//    srand(seed);
//    vector<int> lambda ((unsigned long) 2 * n + 1);
//    lambda[0] = n;
//    for (int i = 1; i <= n; i++) lambda[i] = n;
//    for (int i = n; i <= 2 * n; i++) lambda[i] = 0;
//    vector<int> mu (0);


    BoundedTableaux Top, Bottom, Top_copy, Bottom_copy;

    Top.set_to_top(lambda, 2 * n, t, m);
    Bottom.set_to_bottom(lambda, 2 * n, t, m);
    Top_copy = Top;
    Bottom_copy = Bottom;
    Coupling c(Bottom, Bottom_copy, Top, Top_copy);

    c.cftp(seed, 4096);

    c.result.print_f_vals();
}

void run_skew_shape2() {
    unsigned int seed = 21034;
    int N = 30;
    int t = 30;
    int m = 2;
    int n = N * m;

    vector<int> lambda ((unsigned long) 4 * n + 1);
    lambda[0] = n;
    for (int i = 1; i <= n; i++) lambda[i] = 3*n + 2 - i;
    for (int i = n+1; i <= 2*n; i++) lambda[i] = 2*n + 1;
    for (int i = 2*n+1; i <= 3*n; i++) lambda[i] = 4*n + 1 - i;
    for (int i = 3*n+1; i <= 4*n; i++) lambda[i] = n;
    vector<int> mu (0);


    BoundedTableaux Top, Bottom, Top_copy, Bottom_copy;

    Top.set_to_top(lambda, 4*n, t, m);
    Bottom.set_to_bottom(lambda, 4*n, t, m);
    Top_copy = Top;
    Bottom_copy = Bottom;
    Coupling c(Bottom, Bottom_copy, Top, Top_copy);
//    auto real_rand = bind(std::uniform_real_distribution<double>(0,1),
//                          mt19937(seed));

    c.cftp(seed, pow(2,15));

//    Bottom.print_f_vals();
//    Top.print_f_vals();

//    for (unsigned int i = 0; i < 10; i++) {
//        Bottom.step(real_rand(), rand());
//    }
//    Bottom.print_f_vals();
//    Top.print_f_vals();


    c.result.print_f_vals();


}

void cusps() {
    unsigned int seed = 21034;
    int N = 30;
//    int t = 30;
    int m = 2;
    int n = N * m;
    int t = n;

    vector<int> lambda ((unsigned long) n + 1);
    lambda[0] = n;
    for (int i = 1; i <= n; i++) lambda[i] = (int) floor( i - pow(i,2.0)*1.0/n );
    vector<int> mu (0);


    BoundedTableaux Top, Bottom, Top_copy, Bottom_copy;

    Top.set_to_top(lambda, 5*n, t, m);
    Bottom.set_to_bottom(lambda, 5*n, t, m);
    Top_copy = Top;
    Bottom_copy = Bottom;
    Coupling c(Bottom, Bottom_copy, Top, Top_copy);

    c.cftp(seed, pow(2,15));



    c.result.print_f_vals();


}

void run_skew_shape() {
    unsigned int seed = 21034;
    int N = 60;
    int m = 1;
    int n = N * m;
    int t = n;

    vector<int> lambda ((unsigned long) 2 * n + 1);
    lambda[0] = 8 * n - 1;
    for (int i = 1; i <= n; i++) lambda[i] = 8*n - i;
    for (int i = n+1; i <= 2*n; i++) lambda[i] = 2*n - i + 1;
    vector<int> mu (0);


    BoundedTableaux Top, Bottom, Top_copy, Bottom_copy;

    Top.set_to_top(lambda, 2*n, t, m);
    Bottom.set_to_bottom(lambda, 2*n, t, m);
    Top_copy = Top;
    Bottom_copy = Bottom;
    Coupling c(Bottom, Bottom_copy, Top, Top_copy);
//    auto real_rand = bind(std::uniform_real_distribution<double>(0,1),
//                          mt19937(seed));

    c.cftp(seed, 4096);

//    for (unsigned int i = 0; i < pow(2, 10); i++) {
//        Bottom.step(real_rand(), rand());
//    }
    c.result.print_f_vals();


}

void run_one_path() {
    unsigned int seed = 323144;
    int N = 80;
    int t = 80;
    int m = 2;
    int n = N * m;
    unsigned h = 1;
    int W = n;

    SimpleTableaux Top, Bottom, Top_copy, Bottom_copy;
    Top.W = W; Top.h = h;
    Bottom.W = W; Bottom.h = h;
    Top.set_to_top(n, t, m);
    Bottom.set_to_bottom(n, t, m);
    Top_copy = Top;
    Bottom_copy = Bottom;

    Coupling c(Bottom, Bottom_copy, Top, Top_copy);

    c.cftp(seed, 4096);


    c.result.print_f_vals();

}

void run_nonintersecting_coupling() {
    unsigned int seed = 392322;
    int m = 4;
    int k = 25;
    int j = 25;
    int n = 25;
    int N = k * m;
    int h = j * m / 2;
    int t = 1;

    vector<int> bottom ((unsigned long) N * h);
    vector<int> top ((unsigned long) N * h);

    for (int i = 0; i < h; i++) {
        for (int c = 0; c < N; c++) {
            bottom[i * N + c] = i;
            top[i * N + c] = t * m * n - i;
        }
    }
    vector<int> bottom_copy = bottom, top_copy = top;
    LozengeTiling Bottom, Top, Bottom_copy, Top_copy;
    Bottom.set_params(N, bottom, n, t, m, h);
    Bottom_copy.set_params(N, bottom_copy, n, t, m, h);
    Top.set_params(N, top, n, t, m, h);
    Top_copy.set_params(N, top_copy, n, t, m, h);

    Coupling c(Bottom, Bottom_copy, Top, Top_copy);
    c.cftp(seed, 4096);
    c.result.print();
}

void run_Lozenge_coupling() {
    unsigned int seed = 392322;
    int m = 4;
    int k = 20;
    int j = 20;
    int n = 20;
    int N = k * m;
    int h = j * m / 2;
    int t = 1;

    vector<int> bottom ((unsigned long) N * h);
    vector<int> top ((unsigned long) N * h);

    for (int i = 0; i < h; i++) {
        for (int c = 0; c < N; c++) {
            bottom[i * N + c] = 0;
            top[i * N + c] = t * m * n;
        }
    }
    vector<int> bottom_copy = bottom, top_copy = top;
    LozengeTiling Bottom, Top, Bottom_copy, Top_copy;
    Bottom.set_params(N, bottom, n, t, m, h);
    Bottom_copy.set_params(N, bottom_copy, n, t, m, h);
    Top.set_params(N, top, n, t, m, h);
    Top_copy.set_params(N, top_copy, n, t, m, h);

    Coupling c(Bottom, Bottom_copy, Top, Top_copy);
    c.cftp(seed, 4096);
    c.result.print();
}

void run_bounded() {
    unsigned int seed = 20134;
    int N = 30;
    int t = 30;
    int m = 2;
    int n = N * m;
    double u = 1;
    double v = 1;
    double q = .5;
//    vector<int> tab ((unsigned long) n * n);
    vector<int> lambda ((unsigned long) n+1);
    lambda[0] = 2 * n;
    lambda[1] = 2 * n;
    for (int i = 2; i <= n; i++) lambda[i] = n;
    vector<int> mu (0);


    BoundedTableaux T(lambda, mu, n, u, v, q);

    auto real_rand = bind(std::uniform_real_distribution<double>(0,1),
                          mt19937(seed));

    for (unsigned int i = 0; i < pow(2, 16); i++) {
        T.step(real_rand(), rand());
    }

    T.print_f_vals();

}


void run_simple() {
    unsigned int seed = 26457;
    int t = 20;
    int m = 3;
    int N = 20;
    int n = N * m * 4;
    int W = N * m;
    unsigned long h = (unsigned long) N * m;
    srand(seed);
    vector<int> tab (W * h);

    for (int i = 0; i < h; i++) {
        for (int j = 0; j < W; j++) {
            tab[i * W + j] = t * (n - i) + (t - 1) * j;
        }
    }

    SimpleTableaux T(tab, n, W, h, t, m);

    auto real_rand = bind(std::uniform_real_distribution<double>(0,1),
                          mt19937(seed));

    for (unsigned int i = 0; i < pow(2, 16); i++) {
        T.step(real_rand(), rand());
    }
    // Print out height function restricted to different horizontal lines on the plane
//    cout << '{' << endl;
//
//    for (int i = 2; i < 120; i+=20) {
//        for (unsigned int k = 0; k < n; k++) {
//            int s = 0;
//            for (unsigned int count = 0; count < 200; count++) {
//                s += (T.f(i + 1, k + 1) + T.f(i + 1, k) + T.f(i, k + 1) + T.f(i, k)) / 4.0;
//                T.step(real_rand(), rand());
//            }
//            cout << s / 200.0;
//            if (k < n - 1) {
//                cout << ", " << endl;
//            }
//        }
//        cout << '}' << endl;
//    }
    T.print_f_vals();
}

void run_simple_coupling() {
    unsigned int seed = 26550;
    int t = 2 * 30;
    int m = 2;
    int N = 30;
    int n = N * m;
    int W = N * m;
    unsigned long h = (unsigned long) N * m;


    SimpleTableaux Top, Bottom, Top_copy, Bottom_copy;
    Top.W = W; Top.h = h;
    Bottom.W = W; Bottom.h = h;
    Top.set_to_top(n, t, m);
    Bottom.set_to_bottom(n, t, m);
    Top_copy = Top;
    Bottom_copy = Bottom;

    Coupling c(Bottom, Bottom_copy, Top, Top_copy);

    c.cftp(seed, 4096);
//
//    c.result.print();
//
    c.result.print_f_vals();
}



void print(Tableaux& T) {
    cout << "{" << '\n';
    for (unsigned i=0; i < T.height + 2; ++i) {
        cout << "{";
        for (unsigned j = 0; j < T.N + 2; ++j) {
            cout << T.get(i, j) << ", ";
        }
        cout << "}" << endl;
    }
    cout << "}" << endl;
}

void print_tiling(LozengeTiling& T) {
    cout << "{" << endl;
    for (unsigned i=0; i < T.height + 1; ++i) {
        cout << "{";
        for (unsigned j = 0; j < T.N + 1; ++j) {
            cout << T.get(i, j);
            if (j != T.N) {
                cout << ", ";
            }
        }
        cout << "}";
        if (i != T.height) cout << ",";
        cout << '\n';

    }
    cout << "}" << endl;
}

void print_tiling_derivative(LozengeTiling& T) {
    cout << "{" << endl;
    for (unsigned i=1; i < T.height + 1; ++i) {
        cout << "{";
        for (unsigned j = 1; j < T.N + 1; ++j) {
            if (j != T.N) {
                cout << T.get(i, j) - T.get(i, j + 1);
                cout << ", ";
            } else {
                cout << 0;
            }
        }
        cout << "}";
        if (i != T.height + 1) cout << ",";
        cout << '\n';

    }
    cout << "}" << endl;
}

void run_Lozenge() {
    unsigned int seed = 23943;
    int t = 1;
    int m = 2;
    int n = 25;
    int N = n * m;
    int h = n * m;

    vector<int> tab ((unsigned long) N * h);

    for (int i = 0; i < h; i++) {
        for (int c = 0; c < N; c++) {
            tab[i * N + c] = 0;
        }
    }
    srand(time(0));
    mt19937 mt_rand(time(0));
    auto real_rand = bind(std::uniform_real_distribution<double>(0,1),
                           mt19937(time(0)));

    LozengeTiling T(N, tab, n, t, m, h);
    for (unsigned int i = 0; i < 100000; i++) {
        double r1 = real_rand();
        T.step(r1, rand());
    }
    for (unsigned int i = 0; i < 1000; i++) {
        cout << (T.get(25, 25) + T.get(24, 25) + T.get(25, 24) + T.get(24, 24))/4.0 << ',' << endl;
        double r1 = real_rand();
        T.step(r1, rand());
    }
    T.print();
}

void run_LHT() {
    unsigned int seed = 22134;
    int N = 25;
    int t = 25;
    int m = 2;
    int n = N * m;
    double u = 1;
    double v = 1;
    double q = .5;
    vector<int> tab ((unsigned long) n * n);
    vector<int> lambda ((unsigned long) n);
    for (int i = 0; i < n; i++) lambda.push_back(n);
    vector<int> mu (0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            tab[i * n + j] = t * m * (n + j - i) - j;
        }
    }

    Tableaux T;
    T.set_params(lambda, mu, tab, n, u, v, q);
    print(T);
//    T.step(seed);
    print(T);
}
