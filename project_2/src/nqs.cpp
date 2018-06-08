#include "../include/nqs.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

nqs::nqs() : WaveFunc(){}

void nqs::initialize(){
    double spread = 0.2;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis (0, 1);
    normal_distribution<double> w_dis(0, spread);

    for(int i = 0; i < M; i++){
        a(i) = dis(gen);
    }


    for(int i = 0; i < N; i++){
        b(i) = dis(gen);
    }

    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            W(i,j) = w_dis(gen);
        }

    }
}

double nqs::evaluate(colvec R){

    double exp_term = 0;
    double prod = 1;
    exp_term = exp(accu(-square((R - a))/(2*sigma_2)));
    for(int j = 0; j < N; j++){
    
        double sum_xi_wij = 0;
        for(int l = 0; l < M; l++){
            sum_xi_wij += R(l)*W(l, j);
        }

        prod *= 1 + exp(b(j) + sum_xi_wij/sigma_2);
    }
    return exp_term*prod;
}

double nqs::E_l(colvec R){

    // Calulates the local energy of the given
    // configuration of the system
    //cout << -0.5*laplace(R) << "   " <<  accu(omega_2*square(R))*0.5 << "    " << 0.5*(- laplace(R) + accu(omega_2*square(R))) << endl;  
    return 0.5*(-laplace(R) + omega_2*accu(square(R)));
}

double nqs::E_l_gibbs(colvec R){

    // Calulates the local energy of the given
    // configuration of the system
    return 0.5*(- laplace_gibbs(R) + omega_2*accu(square(R)));
}

double nqs::laplace(colvec R){

    // Calulates the derivatives of the wave function
    // Both first and second derivatives
    // Called upon in local energy function

    double Hj = 0;
    double exp_term = 0;
    double denom = 0;
    double sum_1 = 0;
    double sum_2 = 0;
    double del_ln_psi = 0;
    double del_ln_psi_sq = 0;
    double laplace_psi = 0;
    double sigmoid = 0;
    double sigmoid_deri = 0;
    double laplace_return = 0;
    
    for(int i = 0; i < M; i++){

        for(int j = 0; j < N; j++){

            double sum_xi_wij = 0;
            for(int l = 0; l < M; l++){
                sum_xi_wij += R(l)*W(l, j);
            }
                
            Hj = - b(j) - (sum_xi_wij/sigma_2);
            exp_term = exp(Hj);
            denom = 1 + exp_term;

            sigmoid = 1.0/denom;
            sigmoid_deri = exp_term/(denom*denom);

            sum_1 += W(i,j)*sigmoid;
            sum_2 += W(i,j)*W(i,j)*sigmoid_deri;
        }

        // First and second derivatives of the logarithm
        // of the nqs wave function
        del_ln_psi = -(R(i) - a(i))/sigma_2 + sum_1/sigma_2;
        del_ln_psi_sq = del_ln_psi*del_ln_psi;
        laplace_psi = - 1/sigma_2 + sum_2/sigma_4;

        laplace_return += del_ln_psi_sq + laplace_psi;
    }
    return laplace_return;
}

double nqs::laplace_gibbs(colvec R){

    // Calulates the derivatives of the wave function
    // Both first and second derivatives
    // Called upon in local energy function

    double Hj = 0;
    double exp_term = 0;
    double denom = 0;
    double sum_1 = 0;
    double sum_2 = 0;
    double del_ln_psi = 0;
    double del_ln_psi_sq = 0;
    double laplace_psi = 0;
    double sigmoid = 0;
    double sigmoid_deri = 0;
    double laplace_return = 0;
    double gibbs_factor = 0.5;

    for(int i = 0; i < M; i++){

        for(int j = 0; j < N; j++){

            Hj = - b(j) - sum(R.t()*W.col(j))/sigma_2;
            exp_term = exp(Hj);
            denom = 1 + exp_term;

            sigmoid = 1/denom;
            sigmoid_deri = exp_term/(denom*denom);

            sum_1 += W(i,j)*sigmoid;
            sum_2 += W(i,j)*W(i,j)*sigmoid_deri;
        }

        // First and second derivatives of the logarithm
        // of the nqs wave function
        del_ln_psi = gibbs_factor*(-(R(i) - a(i))/sigma_2 + sum_1/sigma_2);
        del_ln_psi_sq = del_ln_psi*del_ln_psi;
        laplace_psi = gibbs_factor*(-1/sigma_2 + sum_2/sigma_4);

        laplace_return += del_ln_psi_sq + laplace_psi;
    }

    return laplace_return;
}

colvec nqs::drift_force(colvec R){

    // Drift force, F, to be used in importance sampling

    double Hj = 0;
    double exp_j = 0;
    double term = 0;
    double sum_1 = 0;
    double sigmoid = 0;
    colvec grad_psi_sq = colvec(M);
    
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){

            Hj = - b(j) - sum(R.t()*W.col(j))/sigma_2;

            exp_j = exp(Hj);
            term = 1 + exp_j;

            sigmoid = 1/(term);

            sum_1 += W(i,j)*sigmoid;
        }

        // The first derivative of logarithm of the wave function
        grad_psi_sq(i) = -(R(i) - a(i))/sigma_2 + sum_1/sigma_2;
    }
    // F = 2* '(ln(psi))
    colvec drift_force_i = 2*grad_psi_sq;

    return drift_force_i;
}

double nqs::ratio(colvec R, colvec R_p, int k){

    double eval_R = evaluate(R);
    double eval_R_p = evaluate(R_p);
    double prob = (eval_R_p*eval_R_p)/(eval_R*eval_R);
    return prob;
}

void nqs::update_positions(colvec R){
    D = D_p;
}


void nqs::set_params(vec params){

    /*int M_in, int N_in,
    double sigma_in, double sigma_sq_in,
    double sigma_4_in,
    double omega_in, double omega_sq_in,
    double gamma_in
    */
    
    a = colvec(M);
    b = colvec(N);
    
    W = mat(M, N);

    sigma = params[0];
    sigma_2 = params[1];
    sigma_4 = params[2];
    omega = params[3];
    omega_2 = params[4];
    gamma = params[5];

}



