#include <armadillo>
#include "../include/vmc.h"
#include <random>

using namespace arma;
using namespace std;

double NaiveMh::metropolis_hastings(WaveFunc *psi_t){
    mat R_p(size(R));
    R_p = R;
    
    uniform_int_distribution<int> dis_r(0, N_p - 1);
    uniform_real_distribution<double> dis_step(-1, 1);
    uniform_real_distribution<double> dis_p(0, 1);

    int j = dis_r(*gen);
    for(int i = 0; i < N_d; i++){
        R_p(j, i) += dis_step(*gen) * step ;
    }
    
    if(dis_p(*gen) < psi_t -> ratio(R, R_p)){
        R = R_p;
    } 

    return psi_t -> E_l(R);
}

