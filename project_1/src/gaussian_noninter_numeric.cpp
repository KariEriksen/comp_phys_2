#include "../include/gaussian_noninter_numeric.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianNonInterNumeric::GaussianNonInterNumeric() : WaveFunc(){}

double GaussianNonInterNumeric::E_l(mat R){
    double _psi = evaluate(R);
    double _laplace_psi = laplace(R);
    double pot = (double) as_scalar(accu(square(R)));
    
    return - 0.5 * _laplace_psi/_psi + 0.5 * pot ;
}

double GaussianNonInterNumeric::evaluate(mat R){
    double alpha = params[0];
    double beta = params[1];

    if(N_d > 2){
        R.col(2) *= beta;
    }
    double ret_val = 0;
    double internal = accu(square(R));
	
	Dunno why this is here?
    if(N_d >2){
        R.col(2) *= 1/beta;
    }
	
    ret_val = (double) as_scalar(exp(-alpha *(internal)));
    return ret_val;
}

double GaussianNonInterNumeric::laplace(mat R){
    double h = params[2];
    double lap = 0;

    mat Rp = R;
    mat Rm = R;

    for(int i = 0; i<N_p; i++){
        lap += - 2*N_d * evaluate(R);

        for(int j = 0; j < N_d; j++){
           Rp(i, j) += h;
           Rm(i, j) += -h;
            
           double temp = evaluate(Rp) + evaluate(Rm);
           lap += temp;
           
           Rp(i, j) = R(i, j);
           Rm(i, j) = R(i, j);
        }
    } 
    double fact = 1/(h*h);
    return lap*fact;
}

double GaussianNonInterNumeric::drift_force(mat R){
    double h = params[2];
    double der = (evaluate(R + h) - evaluate(R))/h;
    return der;
}


double GaussianNonInterNumeric::ratio(mat R, mat R_p, int k){
    double eval_R = evaluate(R);
    double eval_R_p = evaluate(R_p);

    double prop = eval_R_p / eval_R;
    return prop*prop;
}

void GaussianNonInterNumeric::initialize(mat R){}
void GaussianNonInterNumeric::update(){};

void GaussianNonInterNumeric::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    //params alhpa, beta, h
    params = params_i;
}
