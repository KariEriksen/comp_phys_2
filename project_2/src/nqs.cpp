#include "../include/nqs.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

nqs::nqs() : WaveFunc(){}

void nqs::initialize(mat R, mat a, mat b, mat W){
    uniform_real_distribution<double> dis (0, 0.001);
    gen = new mt19937(rd());

    for(int i = 0; i < N; i++){
        a(i) = dis(*gen);
    }

    for(int j = 0; j < N; j++){
        b(j) = dis(*gen);
    }

    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            W(i,j) = dis(*gen);
        }

    }
}

double nqs::evaluate(mat R){

    double exp_term = 0;
    double prod = 0;
    double sigma_sq = sigma*sigma;
    exp_term = exp(-(R - a)/2*sigma_sq);
    prod *= 1 + exp(b + sum(R*W)/sigma_sq);
    return exp_term*prod;
}

double nqs::E_l(mat R, mat a, mat b, mat W){

    double omega_sq = omega*omega;
    return 0.5*(laplace(R, a, b, W) + sum(omega_sq*R));
}

double nqs::laplace(mat R, mat a, mat b, mat W){

    double Hj = 0;
    double exp_term = 0;
    double term = 0;
    double sum_1 = 0;
    double sum_2 = 0;
    double grad_psi_sq = 0;
    double laplace_psi = 0;
    double sigmoid = 0;
    double sigmoid_deri = 0;

    double sigma_sq = sigma*sigma;
    double sigma_qd = sigma_sq*sigma_sq;

    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            for(int k = 0; k < M; k++){
                Hj = - b(j) - R(i)*W(i,j);
            }
            exp_term = exp(Hj);
            term = 1 + exp_term;

            sigmoid = 1/(term);
            sigmoid_deri = exp_j/(term*term);

            sum_1 += W(i,j)*sigmoid;
            sum_2 += W(i,j)*W(i,j)*sigmoid_deri;
        }

        grad_psi_sq = -(R(i) - a(i))/sigma_sq + sum_1/sigma_sq;
        laplace_psi = -1/sigma_sq + sum_2/sigma_qd;

        laplace_return += (-(grad_psi_sq + laplace_psi) + omega_sq*R(i));
    }

    return laplace_return;
}

mat nqs::drift_force(mat R){

    double Hj = 0;
    double exp_j = 0;
    double term = 0;
    double sum_1 = 0;
    double grad_psi_sq = 0;
    double sigmoid = 0;

    double sigma_sq = sigma*sigma;

    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            for(int k = 0; k < M; k++){
                Hj = - b(j) - R(k)*W(k,j);
            }
            exp_j = exp(Hj);
            term = 1 + exp_j;

            sigmoid = 1/(term);

            sum_1 += W(i,j)*sigmoid;
        }

        grad_psi_sq = -(R(i) - a(i))/sigma_sq + sum_1/sigma_sq;
    }

    mat drift_force_i = 2*deri_psi;

    return drift_force_i;
}

double nqs::ratio(mat R, mat R_p, int k){

    double eval_R = evaluate(R);
    double eval_R_p = evaluate(R_p);
    double prob = (eval_R_p*eval_R_p)/(eval_R*eval_R);
    return prob;
}

void nqs::update_positions(mat R){

    D = D_p;
}

void nqs::update_weights(mat a, mat b, mat W){

    a = a_p;
    b = b_p;
    W = W_p;
}

void nqs::set_params(int M, int N, int N_p, int N_d, double sigma, double omega, double gamma){

    /*
    M_h = M;
    N_x = N;
    N_p = N_p;
    N_d = N_d;
    sigma = sigma;
    omega = omega;
    gamma = gamma;
    */
}



