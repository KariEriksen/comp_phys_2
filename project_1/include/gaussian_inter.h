#pragma once

#include <armadillo>
#include "wavefunc.h"

using namespace arma;
using namespace std; 

class gaussian_inter: public wavefunc{
        double E_l(mat R);
        double evaluate(mat R);
        mat nabla(mat R);
        mat laplace(mat R);
        double ratio(mat R, mat R_p);
};
