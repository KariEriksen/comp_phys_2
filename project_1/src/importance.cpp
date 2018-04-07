#include <armadillo>
#include "../include/vmc.h"

using namespace arma;

double Importance::metropolis_hastings(WaveFunc *psi_t, double prev_E_l){
    mat R_p(size(R));
    R_p = R;
	double dt = 0.05; // dt in [0.001,0.01] should produce stable ground state results.

    uniform_int_distribution<int> dis_r(0, N_p - 1);
    uniform_real_distribution<double> dis_step(-1, 1);
    uniform_real_distribution<double> dis_p(0, 1);
	normal_distribution<double> dis_zeta(0.0, 1.0);
	


    int j = dis_r(*gen);
	double zeta = dis_zeta(*gen);
    mat F_drift(1, N_d);
	F_drift = psi_t -> drift_force(R.row(j));
   
	R_p.row(j) += 0.5*F_drift*dt + zeta*sqrt(dt);

    double P = psi_t -> ratio(R, R_p, j);
    P *= P;

    //R = R_p;
    //psi_t -> update();
    //Do we need to update the wave equation? The drift force
    //is independent of psi, so I think not

	mat F_drift_proposed(1, N_d);
	F_drift_proposed = psi_t -> drift_force(R_p.row(j));

    mat term1(1, N_d);
    mat term2(1, N_d);
    double term3;

    double Green_prev;
    double Green_proposed;

    for(int j = 0; j < N_p; j++){
		term1 = R_p.row(j) - R.row(j) - 0.5*dt*F_drift;
		term2 = R.row(j) - R_p.row(j) - 0.5*dt*F_drift_proposed;
		term3 = 2*dt;

		Green_prev += exp((-(dot(term1,term1)))/term3);
		Green_proposed += exp((-(dot(term2,term2)))/term3);
    }

	double q = (accu(Green_prev)/accu(Green_proposed))*P;
    double eps = dis_p(*gen);
    if(eps < q){
        R = R_p;
        psi_t -> update();
        return psi_t -> E_l(R);
    }
    else{
        return prev_E_l;
    }
}
