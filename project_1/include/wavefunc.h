#pragma once
#include <armadillo>

using namespace std;
using namespace arma;

class WaveFunc{
    public:
        int N_d, N_p; 
        vector<double> params;

        WaveFunc(){}
        virtual void set_params(vector<double> params, int N_d, int N_p) =0;
        virtual void initialize(mat R) = 0;
        virtual void update() = 0;
        virtual ~WaveFunc() {}
    
        virtual double evaluate(mat R) = 0; 
        virtual double E_l(mat R) = 0;
        virtual mat drift_force(mat R) = 0;
        virtual double laplace(mat R) = 0;
        virtual double ratio(mat R, mat R_p, int k) = 0; 
       };
