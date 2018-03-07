#pragma once 

#include <armadillo>
#include "wavefunc.h"

using namespace std;
using namespace arma;

class wavefunc;

class GaussianNonInterNumeric: public wavefunc{
    public: 
        void set_params(vector<double> params);
        double proportion(mat R, mat R_p);
        mat evaluate(mat R);
        mat laplace(mat R);
        mat nabla(mat R);
        //mat NumericalDoubleDerivative(mat R);
        GaussianNonInterNumeric();

};
