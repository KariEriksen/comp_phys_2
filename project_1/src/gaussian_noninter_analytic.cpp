#include "../include/gaussian_noninter_analytic.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianNonInterAnalytic::GaussianNonInterAnalytic() : WaveFunc(){}

double GaussianNonInterAnalytic::E_l(mat R){
    double kin = laplace(R);
    double pot = (double) as_scalar(accu(square(R)));
    return kin + 0.5*pot;
}

double GaussianNonInterAnalytic::evaluate(mat R){
    double alpha = params[0];
	double alpha_sq = params[1];
    double beta = params[2];

    mat R_c(size(R));
    R_c = R;
    if(N_d > 2){
        R_c.col(2) *= beta;
    }
    double ret_val = 0;
    double internal = accu(square(R_c));
    ret_val = (double) as_scalar(exp(-alpha *(internal)));
    return ret_val;
}
double GaussianNonInterAnalytic::laplace(mat R){

    double alpha = params[0];
    double alpha_sq = params[1];
	double beta = params[2];
    double factor = N_d*N_p;
	
	if(N_d == 3){
		R.col(2) = R.col(2)*beta;
	}

    //This is the analytical expression for the second derivative of the
    //wave function. Not calculating the laplacien
    double scnd_der = factor*alpha - 2*alpha_sq*as_scalar(accu(square(R)));
    return scnd_der;
}

vec GaussianNonInterAnalytic::drift_force(mat R){
    return 0;
}


double GaussianNonInterAnalytic::ratio(mat R, mat R_p, int k){
    double eval_R = evaluate(R);
    double eval_R_p = evaluate(R_p);

    double prop =  eval_R_p / eval_R;
    return prop*prop;
}

void GaussianNonInterAnalytic::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    // params : alpha, alpha_squared, beta
    params = params_i;

}

void GaussianNonInterAnalytic::initialize(mat R){
}
void GaussianNonInterAnalytic::update(){
}
