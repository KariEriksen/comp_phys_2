#include "../include/gaussian_noninter_numeric.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianNonInterNumeric::GaussianNonInterNumeric() : WaveFunc(){}

double GaussianNonInterNumeric::E_l(mat R){
    double _psi = eval_g(R);
    double _laplace_psi = laplace(R);
    double pot = (double) as_scalar(accu(sum(square(R))));
    
    return - 0.5 * _laplace_psi/_psi + 0.5 * pot ;
}

double GaussianNonInterNumeric::eval_g(mat R){
    double alpha = params[0];
    double beta = params[1];

    mat R_c(size(R));
    R_c = R;
    if(N_d > 2){
        R_c.col(2) *= beta;
    }
    double ret_val = 0;
    double internal = accu(sum(square(R_c)));
    ret_val = (double) as_scalar(exp(-alpha *(internal)));
    return ret_val;
}
double GaussianNonInterNumeric::laplace(mat R){
    double h = params[2];
    double scnd_der = (eval_g(R-h) - 2* eval_g(R) + eval_g(R + h))/(h*h) ;
    return scnd_der;
}

double GaussianNonInterNumeric::drift_force(mat R){
    double h = params[2];
    double der = (eval_g(R + h) - eval_g(R))/h;
    return der;
}


double GaussianNonInterNumeric::ratio(mat R, mat R_p, int k){
    double eval_R = eval_g(R);
    double eval_R_p = eval_g(R_p);
    return eval_R / eval_R_p;
}

void GaussianNonInterNumeric::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    params = params_i;
}
