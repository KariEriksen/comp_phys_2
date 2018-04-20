#include "../include/gaussian_inter_numeric.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianInterNumeric::GaussianInterNumeric() : WaveFunc(){}

void GaussianInterNumeric::initialize(mat R){
    double a = params[3];

    for(int i = 0; i < N_p; i ++ ){
        mat temp_outer = R.row(i);
        for(int j = (i+1) ; j < N_p; j++){
            mat temp_inner = R.row(j);
            mat temp = temp_inner - temp_outer;
            double dist = (double) sqrt(accu(square(temp)));
            D(i, j) = 1 - a/dist;
        }
    }
}
double GaussianInterNumeric::eval_corr(mat R, int k = -1){
    double a = params[3];
    double ret_val = 1;
    
    if(k != -1){
        D_p = D;

        mat temp_outer = R.row(k);
        for(int i = 0;  k > i ; i++){
            mat temp_inner = R.row(i);
            mat temp = temp_inner - temp_outer;
            double dist = (double) sqrt(accu(square(temp)));
            if(dist > a){
                D_p(i, k) = 1 - a/dist;
            }
            else{
                D_p(i, k) = 0;
            }
        }

        for(int i = k; i < N_p ; i++){
            mat temp_inner = R.row(i);
            mat temp = temp_inner - temp_outer;
            double dist = (double) sqrt(accu(square(temp)));
    
            if(dist > a){
                D_p(k, i) = 1 - a/dist;
            }
            else{
                D_p(k, i) = 0;
            }       
        }

        for(int i = 0; i < N_p; i++){
                for (int j = (i + 1); j < N_p ; j++){
                    ret_val *= D_p(i, j);
                    }
                }
        return ret_val;

    }
    else{
        for(int i = 0; i < N_p; i++){
            for (int j = (i + 1); j < N_p ; j++){
                ret_val *= D(i, j);
                }
            }
    return ret_val;
    }
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
    double tmp = eval_corr(R); 
    return eval_g(R)*tmp;
}

double GaussianInterNumeric::E_l(mat R){

    /*
     * since the sampling guarantees that no state in which any |r_i - rj| 
     * is zero the internal potential is taken to be zero as a consequence.
    '*/
    double tmp = eval_corr(R); 
    double _psi = eval_g(R) *tmp;
    double _laplace_psi = laplace(R);
    double V_ext = 0.5 * (double) as_scalar(accu(sum(square(R))));
    
    return - 0.5 * _laplace_psi/_psi +  V_ext;
}


double GaussianInterNumeric::laplace(mat R){
    double h = params[2];
    double lap = 0;

    mat Rp = R;
    mat Rm = R;

    for(int i = 0; i<N_p; i++){
        for(int j = 0; j < N_d; j++){
           Rp(i, j) += h;
           Rm(i, j) += -h;
            
           double temp = evaluate(Rp) + evaluate(Rm) - 2*evaluate(R) ;
           lap += temp;
           
           Rp(i, j) = R(i, j);
           Rm(i, j) = R(i, j);
        }
    } 
    double fact = 1/(h*h);
    return lap*fact;
}

mat GaussianInterNumeric::drift_force(mat R){
    double h = params[2];
    double der = (evaluate(R + h) - evaluate(R))/h;
    return der;
}


double GaussianInterNumeric::ratio(mat R, mat R_p, int k){
    double prob_R = eval_g(R)*eval_corr(R);
    double prob_R_p = eval_g(R_p)*eval_corr(R_p, k);
    
    double prop = prob_R_p / prob_R;
    return prop*prop;
}

void GaussianInterNumeric::update(){
    D = D_p;
}

void GaussianInterNumeric::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    // alpha, beta, h, a
    params = params_i;
    
    D = mat(N_p, N_p, fill::zeros);
    D_p = mat(N_p, N_p, fill::zeros);
}
