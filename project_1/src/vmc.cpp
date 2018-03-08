#define arma_incl
#ifdef arma_incl
#else
    #include <armadillo>
#endif

#include <random>
#include "../include/vmc.h"

using namespace arma;
using namespace std;

vector<double> vmc::monte_carlo(WaveFunc *psi_t){
    cout << "mc_start" << endl;
    vector<double> E_l; 

    for(int i = 0; i < N_mc; i++){
        E_l.push_back(metropolis_hastings(psi_t));
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

vector<double> vmc::solve(WaveFunc *psi_t){
    // -> is dereferencing and  member access to methods of class.
    cout << "pre ass" << endl;
    uniform_real_distribution<double> dis (0, 1);

    for(int i = 0; i < N_p; i++){
        for(int j = 0; j < N_d; j++){
            //dis and gen must be dereferenced to use
            double tmp;
            tmp = dis(*gen);
            R(i, j) = tmp;
        }
    }
    return monte_carlo(psi_t);
}


