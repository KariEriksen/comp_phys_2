#pragma once
#include <armadillo>

using namespace std;
using namespace arma;

class WaveFunc{
    public:
        int M, N;

        WaveFunc(){}
        virtual void set_params(int M, int N, double gamma) =0;
        virtual void initialize(mat R, mat a_i, mat b_i, mat W_ij) = 0;
        virtual void update(mat a_i, mat b_i, mat W_ij) = 0;
        virtual ~WaveFunc() {}
    
        virtual double evaluate(mat R) = 0; 
        virtual double E_l(mat R) = 0;
        virtual mat drift_force(mat R, int j) = 0;
        virtual double laplace(mat R) = 0;
        virtual double ratio(mat R, mat R_p, int k) = 0; 
       };
