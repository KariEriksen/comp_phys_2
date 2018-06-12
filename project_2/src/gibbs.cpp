#include <armadillo>
#include "../include/vmc.h"
#include <random>

using namespace arma;
using namespace std;

double Gibbs::metropolis_hastings(nqs& psi_t, double prev_E_l){

    vec hj = vec(N);
    double sigma_sq = psi_t.sigma_2;
    uniform_real_distribution<double> dis_h (0,1);
    //Calculate hj from P(xi)

    for(int j = 0; j < N; j++){

        double sum_xi_wij = 0;
        for(int i = 0; i < M; i++){
            sum_xi_wij += R(i)*psi_t.W(i, j);
        }

        double Hj = psi_t.b(j) + (sum_xi_wij/sigma_sq);
        double z = 1/(1 + exp(-Hj));
        /*double z_neg = 1/(1 + exp(Hj));
        double r = dis_h(*gen);
        if(r < z){
            hj(j) = z;
        }
        else{
            hj(j) = z_neg;
        }*/
        hj(j) = dis_h(*gen) < z;
    }

    //Calculate xi+1 from P(hj)
    //accept with probability of one
    for(int i = 0; i < M; i++){

        double sum_hj_wij = 0;
        for(int j = 0; j < N; j++){
            sum_hj_wij += psi_t.W(i, j)*hj(j);
        }
        double my = psi_t.a(i) + sum_hj_wij;
        //Creating a normal dist
        uniform_real_distribution<double> dis(my, sigma_sq);
        R(i) = dis(*gen);
    }

    return psi_t.E_l_gibbs(R);

}


