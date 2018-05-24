#include "../include/nqs.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

nqs::nqs() : WaveFunc(){}

void nqs::initialize(mat R, mat a, mat b, mat W){
    return 0;
}

double nqs::eval_corr(mat R, int k = -1){
    return 0;
}
double nqs::eval_g(mat R){
    return 0;
}

double nqs::evaluate(mat R){
    return 0;
}

double nqs::E_l(mat R){
    return 0;
}

double nqs::laplace(mat R, mat a, mat b, mat W, double omega, double sigma){

    double grad_i_sq = 0;
    double laplace_i = 0;
    double exp_h = 0;
    double omega_sq = omega*omega;
    double sigma_sq = sigma*sigma;

    /*for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            for(int k = 0; k < M; k++){
                sum_1 += R();
            }
            exp_h = - b(j) - sum_1;
        }

        grad_i_sq =

        laplace_return = 0.5*(-(grad_i_sq + laplace_i) + omega_sq*R(i));
    }

    return laplace_return;
    */
    return 0;
}

mat nqs::drift_force(mat R){
    return 0;
}

double nqs::ratio(mat R, mat R_p, int k){
    double eval_R = eval_g(R)*eval_corr(R, k);
    double eval_R_p = eval_g(R_p)*eval_corr(R_p, k);

    double prob = (eval_R_p*eval_R_p)/(eval_R*eval_R);
    return prob;
}

void nqs::update(mat R){
    return 0;
}

void nqs::set_params(int M, int N, double gamma){
    return 0;
}



