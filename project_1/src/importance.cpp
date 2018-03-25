#include <armadillo>
#include "../include/vmc.h"

using namespace arma;

double Importance::metropolis_hastings(WaveFunc *psi_t, double prev_E_l){
    mat R_p(size(R));
    R_p = R;

    uniform_int_distribution<int> dis_r(0, N_p - 1);
    uniform_real_distribution<double> dis_step(-1, 1);
    uniform_real_distribution<double> dis_p(0, 1);

    double F_x = psi_t -> drift_force(R);

    int j = dis_r(*gen);
    for(int i = 0; i < N_d; i++){
        R_p(j, i) += R(j, i) + 0.5*F_x*step + dis_step(*gen)*step;
    }

    double eps = dis_p(*gen);
    double P = psi_t -> ratio(R, R_p, j);
    P *= P;

    //R = R_p;
    //psi_t -> update();
    //Do we need to update the wave equation? The drift force
    //is independent of psi, so I think not
    double F_y = psi_t -> drift_force(R_p);

    double term1 = 0;
    double term2 = 0;
    double term3 = 0;

    double Green_x;
    double Green_y;

    for(int i = 0; i < N_p; i++){
        for(int j = 0; j < N_d; j++){

            term1 = R_p(i, j) - R(i, j) - 0.5*step*F_x;
            term2 = R(i, j) - R_p(i, j) - 0.5*step*F_y;
            term3 = 2*step;

            Green_x += exp((-(term1*term1))/term3);
            Green_y += exp((-(term2*term2))/term3);
        }
    }

    double q = (Green_x/Green_y)*P;

    if(eps < q){
        R = R_p;
        psi_t -> update();
        return psi_t -> E_l(R);
    }
    else{
        return prev_E_l;
    }
}
