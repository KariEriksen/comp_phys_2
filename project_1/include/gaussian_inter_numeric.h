#pragma once

#include <armadillo>
#include "wavefunc.h"

using namespace std;
using namespace arma;

class GaussianInterNumeric: public WaveFunc{
    private: 
        double update_jastrow(mat R);
        double eval_g(mat R);
        double eval_corr(mat R, int k);
        mat D;
        mat D_p;
    public:
        void set_params(vector<double> params, int N_d, int N_p);
        double ratio(mat R, mat R_p, int k );
        double evaluate(mat R);
        double E_l(mat R);
        double laplace(mat R);
        mat drift_force(mat R);
        void update();
        void initialize(mat R);
        double h;
        //mat NumericalDoubleDerivative(mat R);
        GaussianInterNumeric();

};
