#include "../include/gaussian_noninter.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){
    // Derived object of WaveFunc class 
    GaussianNonInterNumeric g;
 
    NaiveMh D;
    double alpha, beta, step, h; 
    int N_p, N_d, N_mc;
    alpha = 1; beta = 1; step = 1; h = 1e-4; 
    N_p = 2; N_d = 2; N_mc = 10;
    
    vector<double> params = {alpha, beta, h};
    g.set_params(params, N_d, N_p);
    //must be called or else you literally have no random numbers
    //args are alpha, beta, N_particles, N_dims, N_mccycles
    D.set_params(alpha, beta, N_p, N_d, N_mc);
    D.step = step;
    //reference must be passed or else you get static linking 
    //
    vector<double> result;
    result = D.solve(&g);

    return 0;
}
