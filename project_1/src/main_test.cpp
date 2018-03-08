#include "../include/gaussian_noninter.h"
#include "../include/vmc.h"

int main(int argc, char *argv[]){
    // Derived object of WaveFunc class 
    GaussianNonInterNumeric g;
    
    naive_mh D;
    double alpha, beta; 
    int N_p, N_d, N_mc;
    alpha = 1; beta = 1;
    N_p = 2; N_d = 2; N_mc = 10;

    //must be called or else you literally have no random numbers
    //args are alpha, beta, N_particles, N_dims, N_mccycles
    D.set_params(alpha, beta, N_p, N_d, N_mc);
    //reference must be passed or else you get static linking 
    D.solve(&g);

    return 0;
}
