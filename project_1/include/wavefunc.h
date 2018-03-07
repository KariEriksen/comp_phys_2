#pragma once
#include <armadillo>

using namespace std;
using namespace arma;

class wavefunc{
    public:
        vector<double> params;

        wavefunc(){}
        virtual void set_params(vector<double> params) =0;

        virtual mat evaluate(mat R) = 0;
        virtual mat nabla(mat R) = 0;
        virtual mat laplace(mat R) = 0;
        virtual double proportion(mat R, mat R_p) = 0; 
       };
