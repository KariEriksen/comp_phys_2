#pragma once

#include <armadillo>
#include "wavefunc.h"

using namespace std;
using namespace arma;

class nqs: public WaveFunc{
private:
    mat D;
    mat D_p;
public:
    int M, N;
    double sigma, sigma_2, sigma_4;
    double omega, omega_2;
    double gamma;
    bool interactive; 

    string id = "nqs";

    colvec a, b;
    mat W;
    colvec drift_force(colvec R);

    void set_params(vec params);
    double ratio(colvec R, colvec R_p, int k);
    double interaction(colvec R, int m, int nd);
    double evaluate(colvec R);
    double E_l(colvec R, int m, int nd);
    double laplace(colvec R);
    double E_l_gibbs(colvec R, int m, int nd);
    double laplace_gibbs(colvec R);
    void initialize();
    void update_positions(colvec R);
    nqs();

};
