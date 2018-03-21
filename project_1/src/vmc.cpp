#define arma_incl
#ifdef arma_incl
#else
    #include <armadillo>
#endif

#include <ctime>
#include <random>
#include "../include/vmc.h"

using namespace arma;
using namespace std;

vector<double> vmc::monte_carlo(WaveFunc *psi_t){
    vector<double> E_l; 
    E_l.push_back(psi_t -> E_l(R));
    for(int i = 1; i < N_mc; i++){
        double tmp = metropolis_hastings(psi_t, E_l[i-1]);
        E_l.push_back(tmp);
    }
   return E_l;
}

void vmc::set_params(double a_in, double b_in, int N, int dim,int mc_cycles){
    a = a_in;
    b = b_in;
    N_p = N;
    N_mc = mc_cycles;
    N_d = dim;
    R = mat(N_p, N_d);
    gen = new mt19937(rd());
}

void vmc::generate_positions(double step_int){
    uniform_real_distribution<double> dis (-1, 1);

    for(int i = 0; i < N_p; i++){
        for(int j = 0; j < N_d; j++){
            //gen must be dereferenced to use
            double tmp = dis(*gen) * step_int;
            R(i, j) = tmp;
        }
    }
 
}

double vmc::solve(WaveFunc *psi_t, string filename){
    // -> is dereferencing and  member access to methods of class.
    
    generate_positions(step);
    psi_t -> initialize(R);
    double evaluated = psi_t -> evaluate(R);
    double step_init = step;
    
    /*
     *Accepting a state in which a particle pair is closer than the permitted
     * distance "a" is unphysical, and a waste of MC-cycles. We guarantee that the
     * first state has a non-zero probability to occur and thus also guarantee that any
     * proposed thate that has a < |r_i - r_j | for any particle pair is rejected
     */
    
    
    while(evaluated == 0){
        step_init += step;
        generate_positions(step_init);
        psi_t -> initialize(R);
        evaluated = psi_t -> evaluate(R);
        cout << step_init << endl;
    }
    
    int start_s = clock();
    vector<double> E_l = monte_carlo(psi_t);
    int end_s = clock();
    double time_spent = (end_s - start_s)/(double (CLOCKS_PER_SEC) * 1000);
    mat arma_e_l = mat(E_l);

    string header = "# N_p: " + to_string(N_p) 
        + "| N_d: " + to_string(N_d) 
        + "| N_mc: " + to_string(N_mc) ;

    ofstream output_file(filename);
    
    output_file << header << endl;
    output_file << time_spent << endl;
    arma_e_l.save(output_file, csv_ascii); 

    output_file.close();
    
    return (double) as_scalar(mean(arma_e_l));
}


