#include "../include/gaussian_noninter_numeric.h"
#include "../include/gaussian_noninter_analytic.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){
    // Derived object of WaveFunc class 
 
    NaiveMh D;
    double beta, step, h; 
    int N_p, N_d, N_mc;
    beta = 1; step = 0.1; h = 1e-5; 
    N_p = 50; N_d = 1; N_mc = pow(2, 20);
    

    for(double alpha = 0.3; alpha < 0.8; alpha += 0.05){
        GaussianNonInterAnalytic g;
        NaiveMh D;
        string sim_type = "NM_NIA";
        vector<double> params = {alpha, beta, h};
        g.set_params(params, N_d, N_p);
        //must be called or else you literally have no random numbers
        //args are alpha, beta, N_particles, N_dims, N_mccycles
        D.set_params(alpha, beta, N_p, N_d, N_mc);
        D.step = step;
        //reference must be passed or else you get static linking 
        //
        double result; 
        string filename = "../data/"+ sim_type+ "_a_" + to_string(alpha) + 
            "_b_" + to_string(beta) +
            "_step_" + to_string(step)+
            ".csv";
        result = D.solve(&g, filename);
        }   
    return 0;
}
