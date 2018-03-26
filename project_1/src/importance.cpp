#include <armadillo>
#include "../include/vmc.h"

using namespace arma;

double Importance::metropolis_hastings(WaveFunc *psi_t, double prev_E_l){
    mat R_p(size(R));
    R_p = R;
	double dt = 0.05 // dt in [0.001,0.01] should produce stable ground state results.

    uniform_int_distribution<int> dis_r(0, N_p - 1);
    uniform_real_distribution<double> dis_step(-1, 1);
    uniform_real_distribution<double> dis_p(0, 1);
	
	normal_distribution<> dis_zeta(0, 1);
	


    int j = dis_r(*gen);
	double zeta = dis_zeta(*gen);
    mat F_drift(1, N_d);
	F_drift = psi_t -> drift_force(R.row(j));
   
	F_drift = psi_t -> drift_force(R(j,i))
	R_p(j) += 0.5*F_drift*dt + zeta*sqrt(dt);

    double eps = dis_p(*gen);
    double P = psi_t -> ratio(R, R_p, j);
    P *= P;

    //R = R_p;
    //psi_t -> update();
    //Do we need to update the wave equation? The drift force
    //is independent of psi, so I think not

	mat F_drift_proposed(1, N_d);
	F_drift_proposed(i) = psi_t -> drift_force(R_p.row(j));

    double term1 = 0;
    double term2 = 0;
    double term3 = 0;

    mat Green_prev(1, N_d);
    mat Green_proposed(1, N_d);

    for(int j = 0; j < N_p; j++){
		term1 = R_p.row(j) - R.row(j, i) - 0.5*dt*F_drift;
		term2 = R.row(j) - R_p.row(j) - 0.5*dt*F_drift_proposed;
		term3 = 2*dt;

		Green_prev += exp((-(term1*term1))/term3);
		Green_proposed += exp((-(term2*term2))/term3);
    }

    mat q(1, N_d) = (Green_prev/Green_proposed)*P;

    if(eps < q){
        R = R_p;
        psi_t -> update();
        return psi_t -> E_l(R);
    }
    else{
        return prev_E_l;
    }
}
