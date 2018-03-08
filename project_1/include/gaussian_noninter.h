#pragma once 

#include <armadillo>
#include "wavefunc.h"

using namespace std;
using namespace arma;

class GaussianNonInterNumeric: public WaveFunc{
    public: 
        void set_params(vector<double> params, int N_d, int N_p);
        double proportion(mat R, mat R_p);
        mat evaluate(mat R);
        mat laplace(mat R);
        mat nabla(mat R);
        //mat NumericalDoubleDerivative(mat R);
        GaussianNonInterNumeric();

};
