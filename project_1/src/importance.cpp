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
	

    mat F_drift(1, N_d);

    int j = dis_r(*gen);
	double zeta = dis_zeta(*gen);
   
	for(int i = 0; i < N_d; i++){
		F_drift(i) = psi_t -> drift_force(R(j,i))
        R_p(j, i) += 0.5*F_x*dt + zeta*sqrt(dt);
    }

    double eps = dis_p(*gen);
    double P = psi_t -> ratio(R, R_p, j);
    P *= P;

    //R = R_p;
    //psi_t -> update();
    //Do we need to update the wave equation? The drift force
    //is independent of psi, so I think not

	mat F_drift_proposed(1, N_d);
	for(int i = 0; i < N_d; i++)
		F_drift_proposed(i) = psi_t -> drift_force(R_p(j,i));

    double term1 = 0;
    double term2 = 0;
    double term3 = 0;

    double Green_x;
    double Green_y;

    for(int i = 0; i < N_p; i++){
		term1 = R_p(j, i) - R(j, i) - 0.5*dt*F_drift(i);
		term2 = R(j, i) - R_p(j, i) - 0.5*dt*F_drift_proposed(i);
		term3 = 2*dt;

		Green_x += exp((-(term1*term1))/term3);
		Green_y += exp((-(term2*term2))/term3);
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
