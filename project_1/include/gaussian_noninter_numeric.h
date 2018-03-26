#pragma once 

#include <armadillo>
#include "wavefunc.h"

using namespace std;
using namespace arma;

class GaussianNonInterNumeric: public WaveFunc{
    public: 
        void set_params(vector<double> params, int N_d, int N_p);
        void update();
        void initialize(mat R);
        double ratio(mat R, mat R_p, int k); 
        double evaluate(mat R);
        double E_l(mat R);
        mat laplace(mat R);
        double drift_force(mat R);
        double h;
        //mat NumericalDoubleDerivative(mat R);
        GaussianNonInterNumeric();

};
