#include <armadillo>
#include "../include/vmc.h"
#include <random>

using namespace arma;
using namespace std;

double NaiveMh::metropolis_hastings(WaveFunc *psi_t){
    cout << "met_h start" << endl;
    mat R_p(size(R));
    R_p = R;

    uniform_int_distribution<int> dis_c(0, N_d - 1);
    uniform_int_distribution<int> dis_r(0, N_p - 1);
    uniform_real_distribution<double> dis_step(0, 1);

    R_p(dis_r(*gen), dis_c(*gen)) += dis_step(*gen) * step ;
    
    if(dis_step(*gen) > psi_t -> proportion(R, R_p)){
        R = R_p;
    } 

    return psi_t -> E_l(R);
}

