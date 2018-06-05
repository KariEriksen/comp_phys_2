#pragma once

#include <armadillo>
#include "wavefunc.h"

using namespace std;
using namespace arma;

class nqs: public WaveFunc{
private:
    mat D;
    mat D_p;
public:
    void set_params(vec params);
    mat a;
    mat b;
    mat W;
    double ratio(mat R, mat R_p, int k);
    double evaluate(mat R);
    double E_l(mat R);
    double laplace(mat R);
    mat drift_force(mat R);
    void initialize();
    void update_positions(mat R);
    void update_weights(mat G);

    double sigma, sigma_2, sigma_4;
    double omega, omega_2;
    double gamma;
    nqs();

};
