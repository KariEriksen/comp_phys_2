#include <armadillo>
#include "../include/vmc.h"

using namespace arma;

double Importance::metropolis_hastings(WaveFunc *psi_t, double prev_E_l){
    mat R_p(size(R));
    R_p = R; 

	// dt in [0.001,0.01] should yield stable values of ground state energy.
	// D from kinetic energy operator <- d_value
	double delta_t = 0.005;
	double d_value = 0.5;
	
	// Which particle to move
	uniform_int_distribution<int> dis_r(0, N_p - 1);
	// How far to move the particle
	uniform_real_distribution<double> dis_step(-1, 1);
	// Epsilon, acceptance parameter
	uniform_real_distribution<double> dis_p(0, 1);

	// Zeta, random gaussian distributed number
	normal_distribution<> dis_zeta(0, 1);
	
	// Pick a random particle
	int j = dis_r(*gen);

	// Calculate drift force for the particle
	mat f_drift(size(R.row(j)));
	f_drift = psi_t -> drift_force(R.row(j));
	
	double zeta = dis_zeta(*gen);
	
	for(int i = 0; i < N_d; i++){
		R_p(j,i) += d_value*f_drift[i]*delta_t + zeta*sqrt(delta_t);
	}

    return psi_t -> ratio(R, R_p, 1);
    
}
