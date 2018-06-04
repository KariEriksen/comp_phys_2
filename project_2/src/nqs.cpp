#include "../include/nqs.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

nqs::nqs() : WaveFunc(){}

void nqs::initialize(mat a, mat b, mat W){
    uniform_real_distribution<double> dis (0, 0.001);
    gen = new mt19937(rd());

    for(int i = 0; i < N; i++){
        a(i) = dis(*gen);
    }

    for(int i = 0; i < N; j++){
        b(i) = dis(*gen);
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

    // Calulates the local energy of the given
    // configuration of the system

    double omega_sq = omega*omega;
    return 0.5*(laplace(R, a, b, W) + sum(omega_sq*R));
}

double nqs::laplace(mat R, mat a, mat b, mat W){

    // Calulates the derivatives of the wave function
    // Both first and second derivatives
    // Called upon in local energy function

    double Hj = 0;
    double exp_term = 0;
    double term = 0;
    double sum_1 = 0;
    double sum_2 = 0;
    double grad_psi = 0;
    double grad_psi_sq = 0;
    double laplace_psi = 0;
    double sigmoid = 0;
    double sigmoid_deri = 0;

    double sigma_sq = sigma*sigma;
    double sigma_qd = sigma_sq*sigma_sq;

    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){

            Hj = - b(j) - sum(R%W.col(j));
            exp_term = exp(Hj);
            term = 1 + exp_term;

            sigmoid = 1/(term);
            sigmoid_deri = exp_j/(term*term);

            sum_1 += W(i,j)*sigmoid;
            sum_2 += W(i,j)*W(i,j)*sigmoid_deri;
        }

        // First and second derivatives of the logarithm
        // of the nqs wave function
        grad_psi = -(R(i) - a(i))/sigma_sq + sum_1/sigma_sq;
        grad_psi_sq = grad_psi*grad_psi;
        laplace_psi = -1/sigma_sq + sum_2/sigma_qd;


        laplace_return += (-(grad_psi_sq + laplace_psi) + omega_sq*R(i));
    }

    return laplace_return;
}

mat nqs::drift_force(mat R, mat a, mat b, mat W){

    // Drift force, F, to be used in importance sampling

    double Hj = 0;
    double exp_j = 0;
    double term = 0;
    double sum_1 = 0;
    double grad_psi_sq = 0;
    double sigmoid = 0;

    double sigma_sq = sigma*sigma;

    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){

            Hj = - b(j) - sum(R%W.col(j));

            exp_j = exp(Hj);
            term = 1 + exp_j;

            sigmoid = 1/(term);

            sum_1 += W(i,j)*sigmoid;
        }

        // The first derivative of logarithm of the wave function
        grad_psi_sq = -(R(i) - a(i))/sigma_sq + sum_1/sigma_sq;
    }

    // F = 2* '(ln(psi))
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

void nqs::update_weights(mat G){

    //Is this sufficient, do we use the same G for all
    //indices?


	// The below code implemented in vmc::gradient_descent
	// but the question above still stands.
	/*
    a += -gamma*G(0);
    b += -gamma*G(1);
    W += -gamma*G(2);
	*/
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



