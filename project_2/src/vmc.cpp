#define arma_incl
#ifdef arma_incl
#else
    #include <armadillo>
#endif

#include <math.h>
#include <stdlib.h>  
#include <ctime>
#include <random>
#include "../include/vmc.h"

using namespace arma;
using namespace std;

void vmc::monte_carlo(WaveFunc *psi_t, metadata *exp_vals){
   exp_vals -> exp_E[0] = psi_t -> E_l(R);
   mat square_R;

   if(compute_extra || compute_obd)
       square_R = square(R);
       if(compute_extra){
            double prod_r_cur = (double) as_scalar(accu(square_R)); 
            exp_vals -> prod_R[0] = prod_r_cur;
            exp_vals -> prod_R_exp_E[0] = prod_r_cur * exp_vals -> exp_E[0];
       }
       if(compute_obd){
            /*r = sqrt(x^2 + y^2)*/
            mat all_radii= sqrt(sum(square_R, 1));
            count_obd(all_radii, exp_vals);
       }

    for(int i = 1; i < N_mc; i++){
        double tmp = metropolis_hastings(psi_t, exp_vals -> exp_E[i-1]);

        if(compute_extra || compute_obd){
            if(exp_vals -> exp_E[i-1] != tmp){
                square_R = square(R);
            } 

            if(compute_extra){
                double prod_r_cur = (double) as_scalar(accu(square_R)); 
                exp_vals -> prod_R[i] = prod_r_cur;
                exp_vals -> prod_R_exp_E[i] = prod_r_cur * tmp;
            }
            if(compute_obd){
                /*r = sqrt(x^2 + y^2)*/
                mat all_radii= sqrt(sum(square_R, 1));
                count_obd(all_radii, exp_vals);
        }
       
        exp_vals -> exp_E[i] = tmp;
    }
    }
}

void vmc::set_params(int N, int dim,int mc_cycles,
        bool meta_bool,
        bool obd_bool
    ){
    N_p = N;
    N_mc = mc_cycles;
    N_d = dim;
    R = mat(N_p, N_d);
    gen = new mt19937(rd());
    compute_extra = meta_bool;
    compute_obd = obd_bool;
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

void vmc::count_obd(mat radii, metadata* all_exp){
    for(int i = 0; i < N_p ; i++){
        double tmp_rad = radii[i];
        for(int j = 0; j < obd_n_bins ; j++){
            double tmp_bin = obd_bins[j];
            if((tmp_rad > (tmp_bin - bin_length)) && (tmp_rad < (tmp_bin + bin_length))){
                all_exp -> obd[j] += 1;
                break;
            }
            if(j == (obd_n_bins-1)) all_exp -> obd[j] ++;
        }
    }

}


vector<double> vmc::solve(WaveFunc *psi_t, string filename){
    // -> is dereferencing and  member access to methods of class.
    
    double outer_limit = 5;
    obd_n_bins = 30*outer_limit;

    generate_positions(step);
    //psi_t -> initialize(R);
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
        //psi_t -> initialize(R);
        evaluated = psi_t -> evaluate(R);
    }
    
    metadata all_exp;
    all_exp.exp_E = new double [N_mc]; 
    
    if(compute_extra){
        all_exp.prod_R = new double [N_mc];
        all_exp.prod_R_exp_E = new double [N_mc];
    }

    if(compute_obd){
        all_exp.obd = new long int[obd_n_bins];
        obd_bins = new double[obd_n_bins];

        for(int i = 0; i<obd_n_bins; i++){
            all_exp.obd[i] = 0;
            obd_bins[i] = 0;
        }

        bin_length = outer_limit / (double)obd_n_bins;
        
        for(int i = 0; i<obd_n_bins; i++){
            obd_bins[i] = (i+1) * bin_length;
        }
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

    
    string header = "#{N_p: " + to_string(N_p) 
        + ", N_d: " + to_string(N_d) 
        + ", N_mc: " + to_string(N_mc)
        + ", alpha: " +to_string(a)
        + ", beta: "+ to_string(b) + "}";

    ofstream output_file("../data/"+filename);
    
    output_file << header << endl;
    arma_e_l.save(output_file, csv_ascii); 
    
    cout << "Writing local energies to file: " << "../data/"+  filename << endl;
    output_file.close();

    if(compute_obd){
        string obd_filename = "../data/obd_"+filename;
        ofstream obd_output(obd_filename);
        double* divisor = new double[obd_n_bins];
        int total = 0;
        
        
        for(int i = 0; i < obd_n_bins; i++){
            total += all_exp.obd[i];
        }
        
        for(int i = 0; i < obd_n_bins ; i++){
            if(N_d == 1){
                divisor[i] = 1;
            }
            else if(N_d == 2){
                divisor[i] = obd_bins[i] ; 
            }
            else{
                divisor[i] = obd_bins[i]*obd_bins[i];
            }
        }


        mat arma_obd = mat(obd_n_bins, 2);
        
        for(int i = 0; i< obd_n_bins; i++){
            arma_obd(i, 1) = all_exp.obd[i]/(double)total; 
            arma_obd(i, 0) = obd_bins[i];
        }
        
        cout << "Writing obd to file: "<< obd_filename << endl;
        obd_output << header << "\n";
        arma_obd.save(obd_output, csv_ascii);
        obd_output.close();
    }

    
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

double vmc::gradient_descent(mat R, mat a, mat b, mat W, double sigma){

    /*double deri_psi_a = 0;
    double deri_psi_b = 0;
    double deri_psi_W = 0;

    double sigma_sq = sigma*sigma;

    for(int i = 0; i < )

    deri_psi_a = (R - a)/sigma_sq;

    for(j = 0; j < N; j++){
        Hj = - b(j) - sum(R%W.col(j));
        deri_psi_b = 1 + exp(Hj);
    }
    */
    return 0;



}
