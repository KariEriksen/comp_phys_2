#include <armadillo>
#include <iostream>
#include <math.h>
#include <random>

using namespace std; 
using namespace arma; 

void to_file(vector<vector<double>> values, double best_E, double best_alpha, int num_dim, int num_particles, int index_best, double best_beta = 1){
    string filename = "../data/vmc_N_"+to_string(num_particles)+"_" + "D_" + to_string(num_dim )+ ".txt";

    cout << "writing to file: " + filename << endl;
    ofstream output_file;
    output_file.open(filename);

    output_file << "# Best_alpha = " + to_string(best_alpha);
    output_file << " Best_beta = " + to_string(best_alpha);
    output_file << " Best_Energy = " + to_string(best_alpha);
    output_file << " index of best_ vals = " + to_string(index_best) << endl;
    output_file << "mean(E_l), std(E_l), Var(E_l), alpha, beta" << endl;
    
    vector<vector<double>>::iterator it; 
    for( it = values.begin(); it != values.end(); it++){
         vector<double> temp = *it ;
        for(int j = 0; j < 5; j++){
            output_file << to_string(temp[j]) + ", ";
        }
        output_file << endl;
    }
    output_file.close();

}

vector<double> output( mat E_l, int N, double alpha){
     mat avg_E_l =  mean(E_l)/N;
     mat std_E_l =  stddev(E_l);
     mat var_E_l =  var(E_l);

    vector<double> ret_val;
    ret_val.push_back(double(as_scalar(avg_E_l)));
    ret_val.push_back(double(as_scalar(std_E_l)));
    ret_val.push_back(double(as_scalar(var_E_l)));
    return ret_val;
    }


double f(double R_i, double R_j){
    return 1;
    }

double local_energy( mat R, double alpha, int N){
    return (- alpha * N * (2* alpha *  accu( sum( square(R))) -  size(R)[1]) +0.5*N* accu( sum( square(R))))  ;
    } 

double prob_dens_ratio( mat R,  mat R_p, double alpha, double beta  = 1){
    if ( size(R)[1] == 3){
             R.col(2) *=  beta;
        }

    return exp(-2*alpha*( accu( sum( square(R))) -  accu( sum( square(R_p)))));
    }

 mat prop_move( mat R, int  N, double step, string sample = "none" ){
    
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    double timestep;
    double D_coeff;

    if(sample == "importance"){
        double sum_r =  accu( sum(R));
        for(int i = 0; i < N; i++){
            R[i] += dis(gen)* sqrt(timestep) * sum_r * timestep * D_coeff ;
            }
        
        }
    else{
        for(int i = 0; i < N; i ++){
             R[i] += step * dis(gen); 
            }
        } 
    return R;
    }

    
vector<double> carlo( mat R,  mat R_p, int n_dim,  int n_carlos, int N, double step,double alpha, double beta = 1) {
    /* Random numbers */
    
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    /*Initialize averages*/
     mat E_l =  mat(n_carlos, 1);

    /* Initialize R*/
     mat R_curr;
    
    for(int i = 0; i < N; i++){
        R[i] = step * dis(gen) ;
        }
    
    /*Initialize the energy and variance to start the monte carlo simulation 
     *  -> choose a trial position R_p  = R + r * step, where r is a random var in [0, 1]
     *  -> Metropolis  algo to accept /reject the step (w  = P(R_p)/P(R))
     *  -> If accepted set R <- R_p
     *  -> Update averages
 
     * */
    
    double ratio, metro_num; 

    for(int i = 0; i < n_carlos; i++){
         R_curr = R;
         R_p = prop_move(R, N, step);
         ratio = prob_dens_ratio(R_curr, R_p, alpha, beta);
         metro_num = dis(gen);

         if(ratio > metro_num){
             R = R_p;
             }

         E_l[i] = local_energy(R, alpha, N);
        }

    return output(E_l, N, alpha);
    }

void variational_mc(double beta_start, double beta_stop, double beta_increment, double alpha_start, double alpha_stop, double alpha_increment,
                    mat R,  mat R_p, int n_dim,  int n_carlos, int N, double step)
{

    double best_alpha = alpha_start; 
    double best_beta = beta_start; 
    double best_E = 1e8;
    int best_index = 0;
    int k = 0;

    vector<double> cur_vals; 
    vector<vector<double> > system_energies; 

    SizeMat size_params = size(R);
    if(size_params[1] == 3){
        for(double beta = beta_start;  beta < beta_stop; beta += beta_increment){
            for(double alpha = alpha_start;  alpha < alpha_increment; alpha += alpha_increment){
                cur_vals = carlo(R, R_p, n_dim, n_carlos, N, step, alpha, beta) ;
                if(cur_vals[0] < best_E){
                    best_E = cur_vals[0];
                    best_alpha = alpha;
                    best_beta = beta; 
                    best_index = k;
                }
                cur_vals.push_back(alpha);
                cur_vals.push_back(beta);

                vector<double> copy_of = cur_vals;
                system_energies.push_back(copy_of);
                k++;
            }
        }
    }
    else{
        for(double alpha = alpha_start;  alpha < alpha_stop; alpha += alpha_increment){
            cur_vals = carlo(R, R_p, n_dim, n_carlos, N, step, alpha) ;
            if(cur_vals[0] < best_E){
                                best_E = cur_vals[0];
                                best_alpha = alpha;
                                best_index = k;
                            }

            cur_vals.push_back(alpha);
            cur_vals.push_back(beta_start);

            vector<double> copy_of = cur_vals;
            system_energies.push_back(copy_of);
            k++;
        }
    }
    to_file(system_energies, best_E, best_alpha, size_params[1], size_params[0], best_index, best_beta);
}

int  main(int argc,  char *argv[])
{
    string output_filename; 
    int n_carlos, N, num_dim;
    int hbar = 1; 
    int e_charge = 1; 
    double alpha, beta; 
    double step;

    if(argc < 3){
        cout << "Wrong usage - run as ./a.out num_cycles num_particles num_dim" << endl  ; 
        exit(1);
        }
    else{
        n_carlos = atoi(argv[1]); 
        N = atoi(argv[2]);
        num_dim = atoi(argv[3]);
     }

    step = 0.1;
    alpha = 0.5;
    beta = 1; 

     mat R =  mat(N, num_dim);
     mat R_p =  mat(N, num_dim);

    //carlo(R, R_p, num_dim, n_carlos, N, step, alpha, beta); 
    variational_mc(0, 1, 0.1, -1, 4, 0.001, R, R_p, num_dim, n_carlos, N, step);
    
    /*
    Initialization steps 
     * Choose an initial vector R with positions for all particles
     * Fix number of monte-carlo steps
     * Choose initial values  for variational parameters and compute |psi_T ^a ( R ) | ^2 
     *
    * */ 
}
