#include <armadillo>
#include "../include/vmc.h"

using namespace arma;

double Importance::metropolis_hastings(WaveFunc *psi_t, double prev_E_l){
    mat R_p(size(R));
    R_p = R;

	double beta = psi_t -> params[2];
	double timestep = 0.05; // timestep in [0.001,0.01] should produce stable ground state results.
	double D = 0.5; // Diffusion coefficient?

	// Set up random number generators
    uniform_int_distribution<int> dis_r(0, N_p - 1);
    uniform_real_distribution<double> dis_step(-1, 1);
    uniform_real_distribution<double> dis_p(0, 1);
	normal_distribution<> dis_zeta(0, 1);
	
    int j = dis_r(*gen);
	double zeta = dis_zeta(*gen);

    mat F_drift(1, N_d);
    mat F_drift_proposed(1, N_d);
	
	// Move one particle j
	R_p.row(j) += 0.5*F_drift*timestep + zeta*sqrt(timestep);
	
	// Scale z-direction with beta
	if(N_d == 3){
		R_p(j,2) *= beta;
	}
	
	// Calculate drift force at current and proposed position
	// ERROR: drift_force returns F_drift as 1x1 matrix for numerical.
	F_drift = psi_t -> drift_force(R.row(j));
	F_drift_proposed = psi_t -> drift_force(R_p.row(j));

	// Calculate the value of greensfunction to be used in metropolis algorithm.
	double GreensFunction = 0.0;

	GreensFunction = dot(0.5*(F_drift_proposed + F_drift),
		(D*timestep*0.5*(F_drift - F_drift_proposed) - R_p.row(j) + R.row(j)));

	// Calculate ratio of proposed and current wavefunction.
    double P = psi_t -> ratio(R, R_p, j);
    P *= P;

	// Calculate factors to compare for metropolis test
	double q = exp(GreensFunction)*P;
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
