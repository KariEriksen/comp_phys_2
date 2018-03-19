#include "../include/gaussian_inter_numeric.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianInterNumeric::GaussianInterNumeric() : WaveFunc(){}

void initialize(mat R){
    double a = params[3];

    D = new trimatu(mat(N_p, N_p).zeros(), 1);
    for(int i = 0; i < N_p; i ++ ){
        mat temp_outer = R.row(i);
        for(int j = 0; j > i; j){
            mat temp_inner = R.row(k);
            D(i, j) = 1 - a/(std::abs((double) sum(temp_outer - temp_inner)));
        }
    }
}

double GaussianInterNumeric::eval_corr(mat R, int k = -1){
    double a = params[3];
    double ret_val = 1;
    D_old = D;
    
    if(k != -1){
        mat r_k = R.row(k);
        for(int i = 0;  k < i ; i++){
            r_ik = std::abs(sum(R.row(i) - r_k));
            if(r_ik > a){
                D(i, k) = r_ik;
            }
            else{
                D(i, k) = 0;
            }
        }
        for(int i = k; i < N_p ; i ){
            r_ik = std::abs(sum(R.row(i) - r_k));
            if(r_ik > a){
                D(k, i) = r_ik;
            }
            else{
                D(k, i) = 0;
            }       
        }
    }

    for(int i = 0; i < N_p; i++){
        for (int j = i + 1; j < N_p -1 ; j++){
            ret_val *= D(i, j);
            }
        }
    }
    return ret_val;
}

double GaussianInterNumeric::E_l(mat R){

    /*
     * since the sampling guarantees that no state in which any |r_i - rj| 
     * is zero the internal potential is taken to be zero as a consequence.
    '*/

    double _psi = eval_g(R);
    double _laplace_psi = laplace(R);
    double V_ext = 0.5 * (double) as_scalar(accu(sum(square(R))));
    double V_int = v_int(R); 
    
    return - 0.5 * _laplace_psi/_psi +  V_ext;
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
double GaussianInterNumeric::laplace(mat R){
    double h = params[2];
    double scnd_der = (eval_g(R-h) - 2* eval_g(R) + eval_g(R + h))/(h*h) ;
    return scnd_der;
}

double GaussianInterNumeric::drift_force(mat R){
    double h = params[2];
    double der = (eval_g(R + h) - eval_g(R))/h;
    return der;
}


double GaussianInterNumeric::ratio(mat R, mat R_p, k){
    double prob_R = eval_g(R)*eval_corr(R);
    double prob_R_p = eval_g(R_p)*eval_corr(R_p, k);
    return prob_R / prob_R_p;
}

void GaussianInterNumeric::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    params = params_i;
}
