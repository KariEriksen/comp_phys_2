#include <armadillo>
#include "../include/vmc.h"

using namespace arma;

double Importance::metropolis_hastings(WaveFunc *psi_t, double prev_E_l){
    mat R_p(size(R));
    R_p = R;
	double dt = psi_t -> params[4];

    uniform_int_distribution<int> dis_r(0, N_p - 1);
    uniform_real_distribution<double> dis_step(-1.0, 1.0);
    uniform_real_distribution<double> dis_p(0.0, 1.0);
	normal_distribution<double> dis_zeta(0.0, 1.0);

	// Pick particle j to move.
    int j = dis_r(*gen);

	// Calculate drift force for that particle.
    mat F_drift(1, N_d);
	F_drift = psi_t -> drift_force(R, j);

	// Move the particle
	double zeta = dis_zeta(*gen);
	R_p.row(j) += 0.5*F_drift*dt + zeta*sqrt(dt);
	
    double P = psi_t -> ratio(R, R_p, j);

	// Calculate drift force for particle j in proposed position
	mat F_drift_proposed(1, N_d);
	F_drift_proposed = psi_t -> drift_force(R, j);


    double Greens = 0.0;

	// Calculate greens functions
	for(int i = 0; i < N_d; i++){
		Greens += 0.5*(F_drift_proposed(i) + F_drift(i))*
			(0.5*dt*0.5*(F_drift(i) - F_drift_proposed(i)) - R_p(j,i) + R(j,i));

	}
	Greens = exp(Greens);

    double eps = dis_p(*gen);
    if(eps < P*Greens){
        R = R_p;
        psi_t -> update();
        return psi_t -> E_l(R);
    }
    else{
        return prev_E_l;
    }
}
