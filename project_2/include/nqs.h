#pragma once

#include <armadillo>
#include "wavefunc.h"

using namespace std;
using namespace arma;

class nqs: public WaveFunc{
private:
    double eval_g(mat R);
    double eval_corr(mat R, int k);
    mat D;
    mat D_p;
public:
    void set_params(int M, int N, double gamma);
    double ratio(mat R, mat R_p, int k);
    double evaluate(mat R);
    double E_l(mat R);
    double laplace(mat R, mat a, mat b, mat W, double omega, double sigma);
    mat drift_force(mat R);
    void initialize(mat R, mat a, mat b, mat W);
    void update(mat R);
    //mat NumericalDoubleDerivative(mat R);
    //??
    nqs();

};
