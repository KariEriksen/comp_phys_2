#include <armadillo>
#include "../include/vmc.h"
#include <random>

using namespace arma;
using namespace std;

double Gibbs::metropolis_hastings(nqs *psi_t, double prev_E_l){

    vec hj = vec(N);
    double sigma_sq = psi_t -> sigma_2;
    //Calculate hj from P(xi)

    for(int j = 0; j < N; j++){

        hj(j) = 1/(1 + exp(- psi_t -> b(j) - sum(R%psi_t -> W.col(j))));
    }

    //Calculate xi+1 from P(hj)
    //accept with probability of one
    for(int i = 0; i < M; i++){

        double my = psi_t -> a(i) + psi_t -> W.row(i)*hj;
        //Creating a normal dist
        uniform_real_distribution<double> dis(my, sigma_sq);
        R(i) = dis(*gen);
    }

    return psi_t -> E_l(R);

}


