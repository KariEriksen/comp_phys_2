#include "../include/vmc.h"
#include "../include/wavefunc.h"
#include "../include/gaussian_noninter_analytic.h"
#include "../include/gaussian_inter_analytic.h"

#include <armadillo>

using namespace std;
using namespace arma;


int main(int argc, char *argv[]){
    double beta, step; 
    double alpha = 0.4;
    int N_p, N_d, N_mc, mc_exp;
 
    if( argc < 3){
        cout << "Wrong usage" << endl;    
        exit(1);
    }
    else{
        N_p = atoi(argv[1]); N_d = atoi(argv[2]); mc_exp = atoi(argv[3]);
    }
    N_mc = pow(2, mc_exp);
    beta = 1;
    step = 0.1;
    
    GaussianNonInterAnalytic g; 
    GaussianInterAnalytic u;

    vector<string> inter_type = {"NIA", "IA"};
    
    for(int i = 0; i < 2; i++){

        NaiveMh D; 

        string filename = "NM_"+inter_type[i]+"_a_"+ 
            to_string(alpha) + 
            "_b_" + to_string(beta) + 
            "_step_" + to_string(step) + 
            "_np_" + to_string(N_p) +
            "_nd_" + to_string(N_d) + ".csv" ;

        vector<double> params = {alpha, alpha*alpha, beta};
        g.set_params(params, N_d, N_p);
        //must be called or else you literally have no random numbers
        //args are alpha, beta, N_particles, N_dims, N_mccycles
        D.set_params(alpha, beta, N_p, N_d, N_mc,
                false, 
                true);
        D.step = step;
        //reference must be passed or else you get static linking 
        //
        vector<double> result; 
        
        if(i == 0) result = D.solve(&g, filename); 
        else result = D.solve(&u, filename);
    }

}

