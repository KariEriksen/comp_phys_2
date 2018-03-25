#include <omp.h>
#include "../include/gaussian_noninter_numeric.h"
#include "../include/gaussian_noninter_analytic.h"
#include "../include/wavefunc.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){
    // Derived object of WaveFunc class 
 
    NaiveMh D;
    double beta, step, h; 
    int N_p, N_d, N_mc, mc_exp;
    beta = 1; step = 0.1; h = 1e-4;
    if( argc < 3){
        cout << "Wrong usage" << endl;    
        exit(1);
    }
    else{
        N_p = atoi(argv[1]); N_d = atoi(argv[2]); mc_exp = atoi(argv[3]);
    }
    N_mc = pow(2, mc_exp);
    
    GaussianNonInterAnalytic g;
    GaussianNonInterNumeric u; 

    vector<vector<double>> numeric_results;
    vector<vector<double>> analytic_results; 
    
    double alpha_start, alpha_end, alpha_step;
    alpha_start = 0.1; alpha_end = 1; alpha_step = 0.05;
    int num_sims = (alpha_end- alpha_start)/alpha_step;
    double *alpha_array = new double[num_sims];
    
    for(int i = 0; i < num_sims ; i++){
        alpha_array[i] = alpha_start + i*alpha_start ;
    }

    for(int i = 0; i < num_sims; i++){
        double alpha = alpha_array[i];
        NaiveMh D;
        string sim_type_a = "NM_NIA";
        vector<double> params = {alpha, alpha*alpha, beta};
        g.set_params(params, N_d, N_p);
        //must be called or else you literally have no random numbers
        //args are alpha, beta, N_particles, N_dims, N_mccycles
        D.set_params(alpha, beta, N_p, N_d, N_mc);
        D.step = step;
        //reference must be passed or else you get static linking 
        //
        vector<double> result; 
        string filename = "../data/"+ sim_type_a+ "_a_" + to_string(alpha) + 
            "_b_" + to_string(beta) +
            "_step_" + to_string(step)+
            "_np_" + to_string(N_p)+
            "_nd_" + to_string(N_d)+
            ".csv";
        result = D.solve(&g, filename);
        analytic_results.push_back(result);
        }
    
     for(int i = 0; i < num_sims; i++){
        double alpha = alpha_array[i];
        NaiveMh E;
        string sim_type_n = "NM_NIN";
        vector<double> params = {alpha, beta, h};
        u.set_params(params, N_d, N_p);
        //must be called or else you literally have no random numbers
        //args are alpha, beta, N_particles, N_dims, N_mccycles
        E.set_params(alpha, beta, N_p, N_d, N_mc);
        E.step = step;
        //reference must be passed or else you get static linking 
        //
        vector<double> result; 
        string filename = "../data/"+ sim_type_n+ "_a_" + to_string(alpha) + 
            "_b_" + to_string(beta) +
            "_step_" + to_string(step)+
            "_np_" + to_string(N_p)+
            "_nd_" + to_string(N_d)+
            ".csv";
        result = E.solve(&u, filename);
        numeric_results.push_back(result);
        }   

    string meta_filename = "../data/NM_NIA_NIN_meta_data_np_" + to_string(N_p)+
        "_nd_" + to_string(N_d)+   
        +".csv";
    ofstream meta_file(meta_filename);
    
    meta_file << "alpha,analytic_energy,analytic_time,numeric_energy,numeric_time" << endl;
    double a = alpha_start;

    for(int i = 0; i<num_sims ; i++){
        meta_file << a << "," << analytic_results[i][0] << "," << analytic_results[i][1] ;
        meta_file << "," << numeric_results[i][0] << "," << numeric_results[i][1] << endl;
        a += alpha_step;
    }

    meta_file.close();
    return 0;
}
