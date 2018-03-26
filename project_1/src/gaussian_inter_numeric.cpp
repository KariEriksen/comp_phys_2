#include "../include/gaussian_inter_numeric.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianInterNumeric::GaussianInterNumeric() : WaveFunc(){}

void GaussianInterNumeric::initialize(mat R){
    double a = params[3];
    D = mat(N_p, N_p);
    D.zeros();
    for(int i = 0; i < N_p; i ++ ){
        mat temp_outer = R.row(i);
        for(int j = (i+1) ; j < N_p; j++){
            mat temp_inner = R.row(j);
            D(i, j) = 1 - a/(std::abs((double) accu(sum(temp_outer - temp_inner))));
        }
    }
}
double GaussianInterNumeric::eval_corr(mat R, int k = -1){
    double a = params[3];
    double ret_val = 1;
    mat D_c = D;

    if(k != -1){
        mat r_k = R.row(k);
        for(int i = 0;  k > i ; i++){
            double r_ik = std::abs(sum(R.row(i) - r_k));
            if(r_ik > a){
                D(i, k) = r_ik;
            }
            else{
                D(i, k) = 0;
            }
        }
        for(int i = k; i < N_p ; i++){
            double r_ik = std::abs(sum(R.row(i) - r_k));
            if(r_ik > a){
                D(k, i) = r_ik;
            }
            else{
                D(k, i) = 0;
            }       
        }
    }
    for(int i = 0; i < N_p; i++){
        for (int j = (i + 1); j < N_p ; j++){
            ret_val *= D(i, j);
            }
        }
    if(k != -1){ 
        D_p = D;
        D = D_c; 
    }
    return ret_val;
}

double GaussianInterNumeric::eval_g(mat R){
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

double GaussianInterNumeric::evaluate(mat R){
    return eval_g(R)*eval_corr(R);
}

double GaussianInterNumeric::E_l(mat R){

    /*
     * since the sampling guarantees that no state in which any |r_i - rj| 
     * is zero the internal potential is taken to be zero as a consequence.
    '*/
    double _psi = evaluate(R);
    double _laplace_psi = laplace(R);
    double V_ext = 0.5 * (double) as_scalar(accu(sum(square(R))));
    
    return - 0.5 * _laplace_psi/_psi +  V_ext;
}


double GaussianInterNumeric::laplace(mat R){
    double h = params[2];
    double scnd_der = (evaluate(R-h)
            - 2* evaluate(R)
            + evaluate(R + h))/(h*h) ;
    return scnd_der;
}

mat GaussianInterNumeric::drift_force(mat R){
    double h = params[2];
	mat der(size(R));
    der = (evaluate(R + h) - evaluate(R))/h;
    return der;
}


double GaussianInterNumeric::ratio(mat R, mat R_p, int k){
    double prob_R = evaluate(R)*eval_corr(R);
    double prob_R_p = evaluate(R_p)*eval_corr(R_p, k);
	
	double prop = prop_R_p / prop_R;
    return prop*prop;
}

void GaussianInterNumeric::update(){
    D = D_p;
}

void GaussianInterNumeric::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    params = params_i;
}
