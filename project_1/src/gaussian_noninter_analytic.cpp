#include "../include/gaussian_noninter_analytic.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianNonInterAnalytic::GaussianNonInterAnalytic() : WaveFunc(){}

double GaussianNonInterAnalytic::E_l(mat R){

    double kin = laplace(R);
    double pot = (double) as_scalar(accu(sum(square(R))));

    return kin + 0.5*pot;
}

double GaussianNonInterAnalytic::evaluate(mat R){
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
double GaussianNonInterAnalytic::laplace(mat R){

    double alpha = params[0];
    double alpha_sq = alpha*alpha;

    double scnd_der = N_d*N_p*alpha - 2*alpha_sq*as_scalar(accu(sum(square(R))));
    return scnd_der;
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
