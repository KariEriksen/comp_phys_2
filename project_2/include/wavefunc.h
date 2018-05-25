#pragma once
#include <armadillo>

using namespace std;
using namespace arma;

class WaveFunc{
public:
    int M, N, N_p, N_d;
    double sigma, omega, gamma;

    WaveFunc(){}
    virtual void set_params(int M, int N, int N_p, int N_d, double sigma, double omega, double gamma) =0;
    virtual void initialize(mat R, mat a, mat b, mat W) = 0;
    virtual void update_positions(mat R) = 0;
    virtual void update_weights(mat a, mat b, mat W) = 0;
    virtual ~WaveFunc() {}
    
    virtual double evaluate(mat R) = 0;
    virtual double E_l(mat R) = 0;
    virtual mat drift_force(mat R, int j) = 0;
    virtual double laplace(mat R, mat a, mat b, mat W) = 0;
    virtual double ratio(mat R, mat R_p, int k) = 0;
};
