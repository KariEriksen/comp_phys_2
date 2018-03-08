#define arma_incl
#ifdef arma_incl
#else
    #include <armadillo>
#endif

#include "../include/vmc.h"

vector<double> vmc::monte_carlo(){
    return vector<double> (1 ,1);
}

void vmc::set_params(double a, double b, int N, int dim,int mc_cycles){
    this -> a = a;
    this -> b = b;
    this -> N_p = N;
    this -> N_mc = mc_cycles;
    this -> N_d = dim;

    this -> R = mat(N_p, N_d);

    mt19937 obj(rd());
    uniform_real_distribution<double> dist(0, 1);
    this -> gen = &obj;
    this -> dis = &dist;
}

void vmc::solve(WaveFunc *psi_t){
    // -> is dereferencing and  member access to methods of class.
    psi_t -> set_params(vector<double> (a, b), N_d, N_p);
    
    cout << "pre ass" << endl;
    for(int i = 0; i < N_p; i++){
        for(int j = 0; j < N_d; j++){
            //dis and gen must be dereferenced to use
            double tmp;
            tmp = (*dis)(gen);
            cout << tmp << endl;

            R(i, j) = tmp;
        }
    }
    
    cout << "c" << endl;
    mat c(N_p, N_d);
    cout << "p_call" << endl;
    c = metropolis_hastings(R, psi_t);
}


