#pragma once

#include <armadillo>
#include "wavefunc.h"

using namespace arma;
using namespace std; 

class gaussian_inter: public wavefunc{
        double E_l(mat R);
        mat nabla(mat R);
        mat laplace(mat R);
        double poportion(mat R, mat R_p);
};
