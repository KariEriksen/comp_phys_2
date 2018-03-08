#include "../include/gaussian_noninter.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianNonInterNumeric::GaussianNonInterNumeric() : WaveFunc(){}

double GaussianNonInterNumeric::E_l(mat R){
    double _over_psi = evaluate(R);
    double _nabla_psi = nabla(R);
    
    return _over_psi*_nabla_psi;
    
}

double GaussianNonInterNumeric::evaluate(mat R){
    double alpha = params[0];
    double beta = params[1];
    
    cout << "alpha beta " <<alpha << " " << beta << endl;
    mat R_c(size(R));
    if(N_d > 2){
        R_c.col(2) *= beta;
    }
    double ret_val = (double) as_scalar(exp(-alpha *(accu(sum(R_c))))); 
    return ret_val; 
}
double GaussianNonInterNumeric::laplace(mat R){
    double h = params[2];
    cout << "h " << h << endl;
    double der = (evaluate(R-h) - 2* evaluate(R) + evaluate(R + h))/(h*h) ; 
    return der;
}

double GaussianNonInterNumeric::nabla(mat R){
    return 0;
}
        

double GaussianNonInterNumeric::proportion(mat R, mat R_p){
    cout << "proportion " << endl ; 
    return evaluate(R) / evaluate (R_p);
}

void GaussianNonInterNumeric::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    params = params_i;
}
