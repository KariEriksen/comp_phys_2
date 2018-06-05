#pragma once
#include <armadillo>

using namespace std;
using namespace arma;

class WaveFunc{
    public:
        int M, N, N_p, N_d;
        mat a, b, W;
        double sigma, omega, gamma;

        WaveFunc(){}
        virtual void set_params(vec params) =0;
        virtual void initialize() = 0;
        virtual void update_positions(mat R) = 0;
        virtual void update_weights(mat G) = 0;
        virtual ~WaveFunc() {}
    
        virtual double evaluate(mat R) = 0;
        virtual double E_l(mat R) = 0;
        virtual mat drift_force(mat R) = 0;
        virtual double laplace(mat R) = 0;
        virtual double ratio(mat R, mat R_p, int k) = 0;
        };
