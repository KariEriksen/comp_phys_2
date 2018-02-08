#include <armadillo>
#include <iostream>
#include <math.h>

using namespace std; 


double g(double  R_i, int n_dim, double alpha, double beta = 1){
    if(n_dim == 1){ 
    return exp(- alpha * R_i) ; 
    }
    }

double f(double R_i, double R_j){
    return 1
    }


arma::mat local_energy(arma::mat R){
    return R;
    } 

arma::mat prob_dens_ratio(arma::mat R, double alpha, double beta  = 1){
    return R;
    }

void carlo(arma::mat R, int n_carlos, int N, double alpha, double beta = 1) {
    
    }

int  main(int argc,  char *argv[])
{
	string output_filename; 
	int n_carlos, N, num_dim;
    int hbar = 1; 
    int e_charge = 1; 
    double alpha, beta;   

    if(argc < 3){
        cout << "Wrong usage - compile as g++ main.cpp num_cycles num_particles" << endl  ; 
        exit(1);
        }
    else{
        n_carlos = atoi(argv[0]); 
        N = atoi(argv[1]);
        num_dim = atoi(argv[2]);
     }

    carlo(R, n_carlos, N, alpha, beta); 

    /*Initialization steps 
     * Choose an initial vector R with positions for all particles
     * Fix number of monte-carlo steps
     * Choose initial values  for variational parameters and compute |psi_T ^a ( R ) | ^2 
     *
     * Initialize the energy and variance to start the monte carlo simulation 
     *  -> choose a trial position R_p  = R + r * step, where r is a random var in [0, 1]
     *  -> Metropolis  algo to accept /reject the step (w  = P(R_p)/P(R))
     *  -> If accepted set R <- R_p
     *  -> Update averages
     * */ 
}
