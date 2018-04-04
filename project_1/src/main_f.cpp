#include "../include/wavefunc.h"
#include "../include/gaussian_noninter_analytic.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){
    NaiveMh D;
    double beta, step, h; 
    int N_p, N_d, N_mc, mc_exp;
    beta = 1; step = 1; h = 1e-4;
    if( argc < 3){
        cout << "Wrong usage" << endl;    
        exit(1);
    }
    else{
        N_p = atoi(argv[1]); N_d = atoi(argv[2]); mc_exp = atoi(argv[3]);
    }
    N_mc = pow(2, mc_exp);
    int N_mc_iter = 1e5; 
    
    GaussianNonInterAnalytic g;

    vector<vector<double>> numeric_results;
    vector<vector<double>> analytic_results; 
    
    double alpha_start = 0.2;
    double gamma;

    int max_sims = 20;
    int sim_n = 0;

    double E_der = 0;
    double alpha = alpha_start;

    double *alpha_array = new double[max_sims];

    while((sim_n < max_sims)){
        NaiveMh D;
        string sim_type_a = "GD_NM_NIA";
        vector<double> params = {alpha, alpha*alpha, beta};
        g.set_params(params, N_d, N_p);
        //must be called or else you literally have no random numbers
        //args are alpha, beta, N_particles, N_dims, N_mccycles
        D.set_params(alpha, beta, N_p, N_d, N_mc_iter, true);
        D.step = step;
        //reference must be passed or else you get static linking 
        //
        vector<double> result;
        string filename = "../data/temp.csv";
        result = D.solve(&g, filename);

        double exp_E = result[0];
        double exp_prod_rsquared = result[1];
        double exp_prod_rsquared_el = result[2];
        double temp_ed = (double) E_der; 
        double tmp_a = (double) alpha; 
        alpha_array[sim_n] = tmp_a;

        E_der = 2*(exp_prod_rsquared_el - exp_E * exp_prod_rsquared);
        if(sim_n == 0) gamma = 1/(100*log10(E_der));

        cout << "-------------------------------------" << endl;
        cout << "exp E: " << exp_E << " | old_alpha: " << alpha << " | gamma: " << gamma<<endl; 

        alpha +=  gamma * E_der;

        if(sim_n>0) gamma = abs((alpha - tmp_a) * (1 / (E_der - temp_ed)));
        
        cout << "new alpha: " << alpha ;
        cout << "| E_der : " << E_der  << endl;
        analytic_results.push_back(result); 
        cout << "-------------------------------------" << endl;

        sim_n ++;
    }

    double min_alpha = 1e8;
    double min_energy = 1e8;

    for(int i = 0; i<max_sims ; i++){
        if(analytic_results[i][0] < min_energy){
            min_alpha = alpha_array[i];
            min_energy = analytic_results[i][0];
        }
    }

    cout << "Optimal alpha: " << min_alpha << endl;
    cout << "Optimal energy: " << min_energy << endl;
    return 0;
}
