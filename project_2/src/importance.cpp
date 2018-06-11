#include <armadillo>
#include "../include/vmc.h"
#include <random>
//
using namespace arma;

double Importance::metropolis_hastings(nqs& psi_t, double prev_E_l){
    colvec R_p(M);
    R_p = R;

    uniform_int_distribution<int> dis_r(0, N_p - 1); // Picks particle 0 or 1
    uniform_real_distribution<double> dis_step(-1, 1);
    uniform_real_distribution<double> dis_p(0.0, 1.0);
    normal_distribution<double> dis_zeta(0, 1);

    int j = dis_r(*gen);
    j = j*N_d;

	// Calculate drift force at all positions.
	
    colvec F_drift(M);
    F_drift = psi_t.drift_force(R);

    // Move only particle j
    double zeta = dis_zeta(*gen);
    
    for (int i = 0; i < N_d; i++){
            R_p(j+i) += 0.5*F_drift(j+i)*dt + zeta*sqrt(dt);
    }

    double P = psi_t.ratio(R, R_p, 1); // 3rd input same as naive_mh.

    // Calculate drift force for configuration in proposed position
    colvec F_drift_proposed(M);
    F_drift_proposed = psi_t.drift_force(R_p);


    double Greens = 0.0;

    // Calculate greens functions
    for (int i = 0; i < M; i++){
	Greens += accu(0.5*(F_drift_proposed(i) + F_drift(i)*
		(0.5*dt*0.5*(F_drift(i) - F_drift_proposed(i)) - R_p(i) + R(i))));
    }
    Greens = exp(Greens);

    double eps = dis_p(*gen);
    if(eps < P*Greens){
        R = R_p;
        //psi_t. update();
        return psi_t.E_l(R);
    }
    else{
        return prev_E_l;
    }
}
