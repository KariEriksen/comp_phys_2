#pragma once
#include <armadillo>

using namespace std;
using namespace arma;

class WaveFunc{
    public:
        string id = "wavefunc";

        WaveFunc(){}
        virtual void set_params(vec params) =0;
        virtual void initialize() = 0;
        virtual void update_positions(colvec R) = 0;
        virtual ~WaveFunc() {}
    
        virtual double evaluate(colvec R) = 0;
        virtual double E_l(colvec R, int m, int nd) = 0;
        virtual colvec drift_force(colvec R) = 0;
        virtual double laplace(colvec R) = 0;
        virtual double E_l_gibbs(colvec R, int m, int nd) = 0;
        virtual double laplace_gibbs(colvec R) = 0;
        virtual double ratio(colvec R, colvec R_p, int k) = 0;
        };
