#define arma_incl
#ifdef arma_incl
#else
    #include <armadillo>
#endif

#include "../include/vmc.h"

vector<double> vmc::monte_carlo(){
    return vector<double> (1 ,1);
}

void vmc::set_params(double a, double b, int N, int mc_cycles){
    this -> a = a;
    this -> b = b;
    this -> N = N;
    this -> mc_cycles = mc_cycles;

    mt19937 obj(rd());
    uniform_real_distribution<double> dist(0, 1);
    this -> gen = &obj;
    this -> dis = &dist;
}

void vmc::solve(WaveFunc *psi_t, bool import){
    // -> is dereferencing and  member access to methods of class.
    psi_t -> params = vector<double> (a, b);
    
    if(import){
        cout << "lick my nuts" << endl;    
    }
    else {
        cout << "or not" << endl;
    };
    
    mat b(2, 2);                                                                                          
    b(0, 0) = 1; b(0, 1) = 2; b(1, 0) = 0;
    mat c(2, 2);
    
    c = metropolis_hastings(b, psi_t);
}


