#ifdef arma_import
    #include <armadillo>
#endif

#include <random>
#include "vmc.h"

using namespace std;
using namespace arma;

vmc::vmc(int n_dim, int N, int n_carlos, double step){
    this ->n_dim = n_dim;
    this ->N = N;
    this ->n_carlos = n_carlos;
    this ->step = step;
}

void vmc::set_simul_params(bool imp, bool grad, bool anal, bool inter){
    this ->importance_sampling = imp;
    this ->gradient_desc = grad;
    this ->anal_derivative = anal; 
    this ->interactive = inter;
}

void vmc::print_params(){
    cout << "-------------------------------------" << endl;
    cout << "Parameters of class instance: " << endl;
    cout << "N_dims: " << n_dim << " || N_particles: " << N << endl;
    cout << "n_carlo cycles: " << n_carlos << " || step: " << step << endl; 
    cout << "importance sampling:  " << importance_sampling << endl;
    cout << "gradient descent:     " << gradient_desc << endl;
    cout << "Analytic derivative : " << anal_derivative << endl;
    cout << "interactive :         " << interactive << endl;
    cout << "-------------------------------------" << endl;
}

void vmc::output_result(){
    string filename = "../data/vmc_N_"+to_string(N)+"_" + "D_" + to_string(n_dim )+ ".txt";

    cout << "writing to file: " + filename << endl;
    ofstream output_file;
    output_file.open(filename);

    output_file << "# Best_alpha = " + to_string(best_alpha);
    output_file << " Best_beta = " + to_string(best_beta);
    output_file << " Best_Energy = " + to_string(best_mean_E);
    output_file << " index of best_ vals = " + to_string(index_best) << endl;
    output_file << "mean(E_l), std(E_l), Var(E_l), alpha, beta" << endl;
    
    vector<vector<double>>::iterator it; 
    for( it = simul_values.begin(); it != simul_values.end(); it++){
         vector<double> temp = *it ;
        for(int j = 0; j < 5; j++){
            output_file << to_string(temp[j]) + ", ";
        }
        output_file << endl;
    }
    output_file.close();
}

int vmc::variational_mc_naive(double beta_start, double beta_stop, double beta_increment,
                double alpha_start, double alpha_stop, double alpha_increment 
)
{
    double best_alpha = alpha_start; 
    double best_beta = beta_start; 
    double best_E = 1e8;
    int best_index = 0;
    int k = 0;

    vector<double> cur_vals; 
    
    if(n_dim == 3){
        for(double beta = beta_start;  beta < beta_stop; beta += beta_increment){
            for(double alpha = alpha_start;  alpha < alpha_increment; alpha += alpha_increment){
                cur_vals = carlo(alpha, beta) ;
                if(cur_vals[0] < best_E){
                    best_mean_E = cur_vals[0];
                    best_alpha = alpha;
                    best_beta = beta; 
                    best_index = k;
                }
                cur_vals.push_back(alpha);
                cur_vals.push_back(beta);

                vector<double> copy_of = cur_vals;
                simul_values.push_back(copy_of);
                k++;
            }
        }
    }
    else{
        for(double alpha = alpha_start;  alpha < alpha_stop; alpha += alpha_increment){
            cur_vals = carlo(alpha, 1) ;
            if(cur_vals[0] < best_E){
                                best_mean_E = cur_vals[0];
                                best_alpha = alpha;
                                best_index = k;
                            }

            cur_vals.push_back(alpha);
            cur_vals.push_back(beta_start);

            vector<double> copy_of = cur_vals;
            simul_values.push_back(copy_of);
            k++;
        }
    } 
}

vector<double> vmc::carlo(double alpha, double beta){
    /* Random numbers */
    
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    /*Initialize averages*/
    mat E_l =  mat(n_carlos, 1);

    /* Initialize R*/
    mat R_curr;
    
    for(int i = 0; i < vmc::n_dim ; i++){
        for(int j = 0; j < N; j++){
            R(i, j) = step * dis(gen) ;
        }
    }
    
    /*Initialize the energy and variance to start the monte carlo simulation 
     *  -> choose a trial position R_p  = R + r * step, where r is a random var in [0, 1]
     *  -> Metropolis  algo to accept /reject the step (w  = P(R_p)/P(R))
     *  -> If accepted set R <- R_p
     *  -> Update averages
 
     * */
    
    double ratio, metro_num; 
    
    if(interactive){
	}
    else{
	for(int i = 0; i < n_carlos; i++){
	     R_curr = R;
	     prop_move();
	     ratio = prob_dens_ratio_noninter(R_curr, alpha, beta);
	     metro_num = dis(gen);

	     if(ratio > metro_num){
		 R = R_p;
		 }

	     E_l[i] = local_energy_noninter(alpha, beta);
	    }
	}
    return output(E_l);
}

void vmc::prop_move(){
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    
    double timestep;
    double D_coeff;

    if(vmc::importance_sampling == true){
        double sum_r =  accu( sum(vmc::R));
        for(int i = 0; i < N; i++){
            R[i] += dis(gen)* sqrt(timestep) * sum_r * timestep * D_coeff ;
            }
        
        }
    else{
        for(int i = 0; i < N; i ++){
            R[i] += step * dis(gen); 
            }
        }
    R_p = R;
}


double vmc::prob_dens_ratio_noninter(mat R_curr, double alpha, double beta){
    if ( n_dim == 3){
            R.col(2) *=  beta;
        }
    return as_scalar(exp(-2*alpha*( accu( sum( square(R_curr))) -  accu( sum( square(R_p))))));
}


double vmc::local_energy_noninter(double alpha, double beta){
    return (- alpha * N * (2* alpha *  accu( sum( square(R))) - 
 		n_dim) +0.5* N* accu( sum( square(R))))  ;
}

vector<double> vmc::output(mat E_l){
    mat avg_E_l =  mean(E_l)/N;
    mat std_E_l =  stddev(E_l);
    mat var_E_l =  var(E_l);

    vector<double> ret_val;
    ret_val.push_back(double(as_scalar(avg_E_l)));
    ret_val.push_back(double(as_scalar(std_E_l)));
    ret_val.push_back(double(as_scalar(var_E_l)));
    return ret_val;
}

int vmc::solve(
        double alpha_start, double alpha_stop, double alpha_increment,
        double beta_start, double beta_stop, double beta_increment){

    vmc::R = mat(N, n_dim);
    vmc::R_p = mat(N, n_dim);
    
    if(interactive){
	        
    }
    
    if(gradient_desc == false)
    {
        variational_mc_naive(alpha_start, alpha_stop, alpha_increment, 
                beta_start, beta_stop, beta_increment);
        return 0;
    }
    else{
        return 1;
    }
    
}



