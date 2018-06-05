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
    void set_params(int M, int N, int N_p, int N_d, double sigma, double omega, double gamma);
    double ratio(mat R, mat R_p, int k);
    double evaluate(mat R);
    double E_l(mat R);
    double laplace(mat R);
    mat drift_force(mat R);
    void initialize();
    void update_positions(mat R);
    void update_weights(mat G);
    //mat NumericalDoubleDerivative(mat R);
    //??
    nqs();

};
