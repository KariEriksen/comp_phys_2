#include "../include/gaussian_inter_numeric.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){
    // Derived object of WaveFunc class 
   double beta, step, h, a; 
    int N_p, N_d, N_mc;
    beta = 1; step = 0.1; h = 1e-5; a = 0.0043; 
    N_p = 10; N_d = 2; N_mc = pow(2, 20);
    
    for(double alpha = 0.3; alpha < 0.8; alpha += 0.05){
        GaussianInterNumeric g;
        NaiveMh D;

        vector<double> params = {alpha, beta, h, a};
        g.set_params(params, N_d, N_p);
        //must be called or else you literally have no random numbers
        //args are alpha, beta, N_particles, N_dims, N_mccycles
        D.set_params(alpha, beta, N_p, N_d, N_mc);
        D.step = step;
        //reference must be passed or else you get static linking 
        //
        double result; 
        string filename = "../data/IN_a_" + to_string(alpha) + 
            "_b_" + to_string(beta) +
            "_step_" + to_string(step)+
            ".csv";
        result = D.solve(&g, filename);
        }
    return 0;
}
