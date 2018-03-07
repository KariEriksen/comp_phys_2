#pragma once

#include <armadillo>
#include "wavefunc.h"

using namespace arma;
using namespace std; 

class gaussian_inter: public wavefunc{
        mat evaluate(mat R);
        mat Derivative(mat R);
        mat DoubleDerivative(mat R);
        double poportion(mat R, mat R_p);
};
