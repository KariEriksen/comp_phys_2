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
    double E_l(mat R, mat a, mat b, mat W);
    double laplace(mat R, mat a, mat b, mat W);
    mat drift_force(mat R);
    void initialize(mat R, mat a, mat b, mat W);
    void update(mat R, mat a, mat b, mat W);
    //mat NumericalDoubleDerivative(mat R);
    //??
    nqs();

};
