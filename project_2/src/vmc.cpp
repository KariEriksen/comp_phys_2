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

void vmc::monte_carlo(nqs *psi_t, metadata *exp_vals){
    exp_vals -> exp_E[0] = psi_t -> E_l(R);
    
    mat square_R;

    if(compute_obd)
        square_R = square(R);
    if(compute_extra){
        gradient_descent(psi_t, exp_vals, exp_vals -> exp_E[0]);
    }
    if(compute_obd){
        /*r = sqrt(x^2 + y^2)*/
        mat all_radii= sqrt(sum(square_R, 1));
        count_obd(all_radii, exp_vals);
    }

    for(int i = 1; i < N_mc; i++){
        double tmp = metropolis_hastings(psi_t, exp_vals -> exp_E[i-1]);
        exp_vals -> exp_E[i] = tmp;
        
        if(compute_extra || compute_obd){
            if((exp_vals -> exp_E[i-2] != tmp) & compute_obd){
                square_R = square(R);
            }

            if(compute_extra){
                gradient_descent(psi_t, exp_vals, tmp);
            }
            if(compute_obd){
                /*r = sqrt(x^2 + y^2)*/
                mat all_radii= sqrt(sum(square_R, 1));
                count_obd(all_radii, exp_vals);
            }

        }
    }
}

void vmc::set_params(int N_in, int M_in, int mc_cycles,
                     bool meta_bool,
                     bool obd_bool
                     ){

    gradient_a = mat(M_in, 1);
    gradient_b = mat(N_in, 1);
    gradient_w = mat(M_in, N_in);
    
    N_mc = mc_cycles;
    R = colvec(M_in);

    M = M_in;
    N = N_in;
    gen = new mt19937(rd());
    compute_extra = meta_bool;
    compute_obd = obd_bool;
}

void vmc::generate_positions(double step_int){
    uniform_real_distribution<double> dis (-1, 1);

    for(int i = 0; i < M; i++){
        //gen must be dereferenced to use
        double tmp = dis(*gen) * step_int;
        R(i) = tmp;
    }
}

void vmc::count_obd(mat radii, metadata* all_exp){
    for(int i = 0; i < M ; i++){
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


retval vmc::solve(nqs *psi_t, string filename){
    // -> is dereferencing and  member access to methods of class.
    
    double outer_limit = 5;
    obd_n_bins = 30*outer_limit;

    generate_positions(step);

    /*
     *Accepting a state in which a particle pair is closer than the permitted
     * distance "a" is unphysical, and a waste of MC-cycles. We guarantee that the
     * first state has a non-zero probability to occur and thus also guarantee that any
     * proposed thate that has a < |r_i - r_j | for any particle pair is rejected
     */
    
   
    metadata all_exp;
    all_exp.exp_E = new double[N_mc];
    
    if(compute_extra){
        all_exp.grad_a = colvec(M);
        all_exp.grad_b = colvec(N);
        all_exp.grad_W = mat(M, N);

        all_exp.prod_E_grad_a = colvec(M);
        all_exp.prod_E_grad_b = colvec(N);
        all_exp.prod_E_grad_W = mat(M, N);
        
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
    
    double omc = 1/(double) N_mc;

    all_exp.grad_a *= omc;
    all_exp.prod_E_grad_a *= omc;

    all_exp.grad_b *= omc;
    all_exp.prod_E_grad_b *= omc;

    all_exp.grad_W *= omc;
    all_exp.prod_E_grad_W *= omc;

    mat arma_e_l = mat(N_mc, 1);

    for(int i = 0; i < N_mc; i++) arma_e_l(i, 0) = all_exp.exp_E[i];
    
    /*string header = "#{N_p: " + to_string(N_p)
            + ", N_d: " + to_string(N_d)
            + ", N_mc: " + to_string(N_mc)
            + ", alpha: " +to_string(a)
            + ", beta: "+ to_string(b) + "}";

    ofstream output_file("../data/"+filename);
    
    output_file << header << endl;
    arma_e_l.save(output_file, csv_ascii);
    
    cout << "Writing local energies to file: " << "../data/"+  filename << endl;
    output_file.close();
    */

    if(compute_obd){
        string obd_filename = "../data/obd_"+filename;
        ofstream obd_output(obd_filename);
        double* divisor = new double[obd_n_bins];
        int total = 0;
        
        
        for(int i = 0; i < obd_n_bins; i++){
            total += all_exp.obd[i];
        }

        int N_d = M/2 ;
        
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
        
        /*cout << "Writing obd to file: "<< obd_filename << endl;
        obd_output << header << "\n";
        arma_obd.save(obd_output, csv_ascii);
        obd_output.close();
        */
    }

    
    retval r;
    r.exp_vals = all_exp;
    r.time_spent = time_spent;
    r.el_exp = accu(arma_e_l)/(double) N_mc ; 
    return r;
}

void vmc::gradient_descent(nqs *psi_t, metadata *exp_vals, double E_l){

    // Visible biases are vectors of length M.
    // Hidden nodes and corresponding hidden biases are vectors of length N
    // Weight matrix MxN
    // Number of visible nodes should correspond to the number of particles, N_p,
    // and the number of dimensions in the system, N_d, such that
    // M = N_p*N_d

    //int M = R.size(0);
    //int N = W.size(1);
    

   double sigma_squared = psi_t -> sigma_2;


    // Gradient a
    gradient_a += R - psi_t -> a;
    gradient_a /= sigma_squared;

    // Gradient b
    for (int k = 0; k < N; k++) {
        double inner_sum;
        inner_sum = accu(R.t()*psi_t -> W.col(k))/sigma_squared;
        gradient_b(k) = 1/(1 + exp(-psi_t -> b(k) - inner_sum));
    }

    // Gradient w_kn
    // Concider changing summation indices in .tex for clarity
    for (int n = 0; n < N; n++) {
        for (int k = 0; k < M; k++) {
            double temp = 0.0;
            temp += accu(R.t()*psi_t->W.col(n))/sigma_squared;
            gradient_w(k, n) += 1/(1 + exp(-psi_t->b(n) - temp));
        }
    }
    gradient_w /= sigma_squared;

    exp_vals -> grad_a += gradient_a;
    exp_vals -> grad_b += gradient_b;
    exp_vals -> grad_W += gradient_w;

    exp_vals -> prod_E_grad_a += E_l * gradient_a;
    exp_vals -> prod_E_grad_b += E_l * gradient_b;
    exp_vals -> prod_E_grad_W += E_l * gradient_w;
}
