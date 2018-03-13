#include "../include/gaussian_noninter_analytic.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianNonInterAnalytic::GaussianNonInterAnalytic() : WaveFunc(){}

double GaussianNonInterAnalytic::E_l(mat R){
    double _over_psi = evaluate(R);
    double _laplace_psi = laplace(R);

    return _over_psi*_laplace_psi;

}

double GaussianNonInterAnalytic::evaluate(mat R){
    double alpha = params[0];
    double beta = params[1];

    mat R_c(size(R));
    R_c = R;
    if(N_d > 2){
        R_c.col(2) *= beta;
    }
    double ret_val = (double) as_scalar(N_d*N_p - 2*alpha*alpha*(sum(R_c)));
    //double ret_val = (double) as_scalar(exp(-alpha *(accu(sum(R_c)))));
    return ret_val;
}
double GaussianNonInterAnalytic::laplace(mat R){
    return 0;
}

double GaussianNonInterAnalytic::nabla(mat R){
    return 0;
}


double GaussianNonInterAnalytic::ratio(mat R, mat R_p){
    double eval_R = evaluate(R);
    double eval_R_p = evaluate(R_p);

    double prop = eval_R / eval_R_p;
    return prop;
}

void GaussianNonInterAnalytic::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    params = params_i;
}
