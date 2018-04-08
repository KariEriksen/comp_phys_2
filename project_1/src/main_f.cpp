#include "../include/wavefunc.h"
#include "../include/gaussian_noninter_analytic.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){
    NaiveMh D;
    double beta, step, h; 
    int N_p, N_d, N_mc, mc_exp;
    double alpha_start;

    beta = 1; step = 1; h = 1e-4;
    if( argc < 3){
        cout << "Wrong usage" << endl;    
        exit(1);
    }
    else{
        N_p = atoi(argv[1]); N_d = atoi(argv[2]); mc_exp = atoi(argv[3]);
        alpha_start = atof(argv[4]);
    }
    N_mc = pow(2, mc_exp);
    int N_mc_iter = 1e5; 
    
    GaussianNonInterAnalytic g;

    vector<vector<double>> numeric_results;
    vector<vector<double>> analytic_results; 
    
    double gamma = 1;
    double s_km, y_km;

    int max_sims = 40;
    int k = 0;

    double E_der = 0;
    double epsilon = 1e-6;

    double *alpha_array = new double[max_sims];
    double *gamma_array = new double[max_sims];
    double *gradient_array = new double[max_sims];

    alpha_array[0] = alpha_start;
    
    while(k < max_sims){
        NaiveMh D;
        string sim_type_a = "GD_NM_NIA";
        vector<double> params = {alpha_array[k], alpha_array[k]*alpha_array[k], beta};
        g.set_params(params, N_d, N_p);
        //must be called or else you literally have no random numbers
        //args are alpha, beta, N_particles, N_dims, N_mccycles
        D.set_params(alpha_array[k], beta, N_p, N_d, N_mc_iter, true);
        D.step = step;
        //reference must be passed or else you get static linking 
        vector<double> result;
        string filename = "../data/temp.csv";
        result = D.solve(&g, filename);
        double exp_E = result[0];
        double exp_prod_rsquared = result[1];
        double exp_prod_rsquared_el = result[2];

        cout << "-------------------------" << endl;
        E_der = 2*(exp_prod_rsquared_el - exp_E * exp_prod_rsquared);
 
        gradient_array[k] = (double) E_der;

        if(abs(E_der) < epsilon) break;

        if(k == 0) gamma_array[k] = E_der/1000;
        else{
            s_km = alpha_array[k] - alpha_array[k-1];
            y_km = gradient_array[k] - gradient_array[k-1];
            gamma_array[k] = s_km/y_km;
        }  
        
        cout << "a " << alpha_array[k] << " | g " << gradient_array[k] << endl; 
        cout << "da " << s_km << endl;  
        cout << "de " << y_km << endl;
        cout << "gamma " << gamma_array[k] << endl;
        cout << "iter " << k << endl;

        alpha_array[k+1] = alpha_array[k] - gamma_array[k] * E_der;
        if(abs(alpha_array[k+1] - alpha_array[k]) < epsilon) break;

        if(alpha_array[k+1] < 1e-2){
            alpha_array[k+1] = 1e-1;
            gamma_array[k] = (alpha_array[k] - alpha_array[k+1])/E_der;
        }

        analytic_results.push_back(result); 
        
        k ++;
    }
    double min_alpha = 1e8;
    double min_energy = 1e8;
    
    int j = 0;
    for(auto it = begin(analytic_results); it != end(analytic_results); ++it){
        if(it[0][0] < min_energy) {
            min_alpha = alpha_array[j];
            min_energy = it[0][0];
        }
        j++;
    }
    
    string filename = "../data/GD_NM_NIA_np_"+to_string(N_d)+"_nd_"+to_string(N_d);
    filename += "_sa_"+to_string(alpha_start)+".csv"; 
    ofstream fo(filename);
    for(int i = 0; i<j; i++)fo << alpha_array[i] << "\n";
    fo.close();

    return 0;
}
