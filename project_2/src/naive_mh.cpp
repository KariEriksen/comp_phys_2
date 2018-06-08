#include <armadillo>
#include "../include/vmc.h"
#include <random>

using namespace arma;
using namespace std;

double NaiveMh::metropolis_hastings(nqs *psi_t, double prev_E_l){
    colvec R_p(M);
    R_p = R;

    uniform_int_distribution<int> dis_r(0, N_p-1);
    uniform_real_distribution<double> dis_step(-1, 1);
    uniform_real_distribution<double> dis_p(0, 1);

    int j = dis_r(*gen);
    for(int i = 0; i < N_p * N_d; i++){
        R_p(j*i) += dis_step(*gen) * step ;
    }

    double eps = dis_p(*gen);
    double P = psi_t -> ratio(R, R_p, 1);

    if(eps < P){
        R = R_p;
        psi_t -> update_positions(R);
        return psi_t -> E_l(R);
    }
    else{
        return prev_E_l;
    }

}


