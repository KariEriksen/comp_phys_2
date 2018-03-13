#pragma once 

#include <armadillo>
#include "wavefunc.h"

using namespace std;
using namespace arma;

class GaussianNonInterNumeric: public WaveFunc{
    public: 
        void set_params(vector<double> params, int N_d, int N_p);
        double ratio(mat R, mat R_p); 
        double evaluate(mat R);
        double E_l(mat R);
        double laplace(mat R);
        double nabla(mat R);
        double h;
        //mat NumericalDoubleDerivative(mat R);
        GaussianNonInterNumeric();

};
