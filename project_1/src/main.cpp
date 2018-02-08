#include <armadillo>
#include <iostream>
#include <math.h>
#include <random>

using namespace std; 


double g(double  R_i, int n_dim, double alpha, double beta = 1){
    if(n_dim == 1){ 
        return exp(- alpha * R_i) ; 
    }
    }

double f(double R_i, double R_j){
    return 1;
    }


arma::mat local_energy(arma::mat R){
    return R;
    } 

arma::mat prob_dens_ratio(arma::mat R, arma::mat R_p, double alpha, double beta  = 1){
    return exp(-2*alpha*(arma::sum(R) - arma::sum(R_p)));
    }

arma::mat prop_move(arma::mat R, int  N, double step){
    
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    for(int i = 0; i < N; i ++){
         R[i] += step * dis(gen); 
        } 
    return R;
    }

void carlo(arma::mat R, arma::mat R_p, int n_dim,  int n_carlos, int N, double step,double alpha, double beta = 1) {
    cout << step << endl << "this was step" << endl; 
    /* Random numbers */
    
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    /*Initialize averages*/
    double E_l = 0;

    /* Initialize R*/
    arma::mat R_curr;
    
    for(int i = 0; i < N; i++){
        R[i] = step * dis(gen) ;
        }
    
    /*Initialize the energy and variance to start the monte carlo simulation 
     *  -> choose a trial position R_p  = R + r * step, where r is a random var in [0, 1]
     *  -> Metropolis  algo to accept /reject the step (w  = P(R_p)/P(R))
     *  -> If accepted set R <- R_p
     *  -> Update averages
 
     * */
    
    for(int i = 0; i < n_carlos; i++){
         R_curr = R;
         R_p = prop_move(R, N, step);
        }
    }

int  main(int argc,  char *argv[])
{
	string output_filename; 
	int n_carlos, N, num_dim;
    int hbar = 1; 
    int e_charge = 1; 
    double alpha, beta; 
    double step;

    if(argc < 3){
        cout << "Wrong usage - compile as g++ main.cpp num_cycles num_particles" << endl  ; 
        exit(1);
        }
    else{
        n_carlos = atoi(argv[1]); 
        N = atoi(argv[2]);
        num_dim = atoi(argv[3]);
     }
    
    step = 0.1;
    alpha = 1;
    beta = 1; 

    arma::mat R = arma::mat(N, num_dim);
    arma::mat R_p = arma::mat(N, num_dim);

    carlo(R, R_p, num_dim, n_carlos, N, step, alpha, beta); 

    /*Initialization steps 
     * Choose an initial vector R with positions for all particles
     * Fix number of monte-carlo steps
     * Choose initial values  for variational parameters and compute |psi_T ^a ( R ) | ^2 
     *
    * */ 
}
