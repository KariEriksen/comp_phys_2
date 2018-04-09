#define arma_incl
#ifdef arma_incl
#else
    #include <armadillo>
#endif

#include <ctime>
#include <random>
#include "../include/vmc.h"

using namespace arma;
using namespace std;

void vmc::monte_carlo(WaveFunc *psi_t, metadata *exp_vals){
   exp_vals -> exp_E[0] = psi_t -> E_l(R);
   int c = 0;

    for(int i = 1; i < N_mc; i++){
        double tmp = metropolis_hastings(psi_t, exp_vals -> exp_E[i-1]);
        exp_vals -> exp_E[i] = tmp;
        if(tmp == exp_vals -> exp_E[i-1]){
            c ++;
        }
        if(compute_extra){
            double prod_r_cur = (double) as_scalar(accu(square(R)));
            exp_vals -> prod_R[i] = prod_r_cur;
            exp_vals -> prod_R_exp_E[i] = prod_r_cur * tmp;
        }
    }
    cout << "prosent" << "   " << (double) c/N_mc << endl;
}

void vmc::set_params(double a_in, double b_in, 
        int N, int dim,int mc_cycles,
        bool meta_bool
    ){
    a = a_in;
    b = b_in;
    N_p = N;
    N_mc = mc_cycles;
    N_d = dim;
    R = mat(N_p, N_d);
    gen = new mt19937(rd());
    compute_extra = meta_bool;
}

void vmc::generate_positions(double step_int){
    uniform_real_distribution<double> dis (-1, 1);

    for(int i = 0; i < N_p; i++){
        for(int j = 0; j < N_d; j++){
            //gen must be dereferenced to use
            double tmp = dis(*gen) * step_int;
            R(i, j) = tmp;
        }
    }
 
}

vector<double> vmc::solve(WaveFunc *psi_t, string filename){
    // -> is dereferencing and  member access to methods of class.
    
    generate_positions(step);
    psi_t -> initialize(R);
    double evaluated = psi_t -> evaluate(R);
    double step_init = step;
    
    /*
     *Accepting a state in which a particle pair is closer than the permitted
     * distance "a" is unphysical, and a waste of MC-cycles. We guarantee that the
     * first state has a non-zero probability to occur and thus also guarantee that any
     * proposed thate that has a < |r_i - r_j | for any particle pair is rejected
     */
    
    while(evaluated == 0){
        step_init += step*1e-5;
        generate_positions(step_init);
        psi_t -> initialize(R);
        evaluated = psi_t -> evaluate(R);
    }
    
    metadata all_exp;
    all_exp.exp_E = new double [N_mc]; 
    
    if(compute_extra){
        all_exp.prod_R = new double [N_mc];
        all_exp.prod_R_exp_E = new double [N_mc];
    }

    int start_s = clock();
    monte_carlo(psi_t, &all_exp);
    int end_s = clock();

    double time_spent = (end_s - start_s)/(double (CLOCKS_PER_SEC) * 1000);

    mat arma_e_l = mat(N_mc, 1);
    mat arma_prod_r= mat(N_mc, 1);
    mat arma_prod_r_el= mat(N_mc, 1);

    for(int i = 0; i < N_mc; i++) arma_e_l(i, 0) = all_exp.exp_E[i];
    
    if(compute_extra){
        for(int i = 0; i < N_mc; i++) arma_prod_r(i, 0) = all_exp.prod_R[i];
        for(int i = 0; i < N_mc; i++) arma_prod_r_el(i, 0) = all_exp.prod_R_exp_E[i];
    } 

    string header = "# N_p: " + to_string(N_p) 
        + "| N_d: " + to_string(N_d) 
        + "| N_mc: " + to_string(N_mc) ;

    ofstream output_file(filename);
    
    output_file << header << endl;
    arma_e_l.save(output_file, csv_ascii); 

    output_file.close();
    
    vector<double> retval;
    if(compute_extra){
         retval = {(double) as_scalar(mean(arma_e_l)),
                                  (double) as_scalar(mean(arma_prod_r)),
                                  (double) as_scalar(mean(arma_prod_r_el)),
                                  time_spent
                                 };
    }
    else{
        retval = {(double) as_scalar(mean(arma_e_l)), time_spent};
    }
    return retval ;
}

