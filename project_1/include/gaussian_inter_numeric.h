#pragma once

#include <armadillo>
#include "wavefunc.h"

using namespace std;

class GaussianInterNumeric: public WaveFunc{
    private: 
        double update_jastrow(mat R);
        double eval_g(arma::mat R);
        double eval_corr(arma::mat R, int k);
        arma::mat D;
        arma::mat D_p;

    public:
        void set_params(vector<double> params, int N_d, int N_p);
        double ratio(arma::mat R, arma::mat R_p, int k );
        double evaluate(arma::mat R);
        double E_l(arma::mat R);
        double laplace(arma::mat R);
        double drift_force(arma::mat R);
        void update();
        void initialize(arma::mat R);
        double h;
        //mat NumericalDoubleDerivative(mat R);
        GaussianInterNumeric();

};
