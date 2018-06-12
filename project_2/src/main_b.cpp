#include "../include/nqs.h"
#include "../include/wavefunc.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

    int N, M;
    int N_p, N_d, N_mc, mc_exp;
	double omg, sigma, step;


    N_p = 2;
    N_d = 3;
    mc_exp = 17;
    N = 2;

    N_mc = pow(2, mc_exp);
    M = N_p*N_d;

    //gamma = 1e-2; 
    double gamma_vals[] = {1e-3, 1e-2, 1e-1, 0.4};
    ofstream timefile("../data/b_data/time_iter.csv");

	// Initializations:
	// We store initialization values for positions, biases, and weights,
	// to improve comparison between different gamma values.
	omg = 1; sigma = 1; step = 0.8;
	double omg_2 = omg*omg;
	double sigm_2 = sigma*sigma;
	double sigm_4 = sigm_2*sigm_2;
	
	NaiveMh D;
	D.step = step;
	D.set_params(N_p, N_d, N, M, N_mc, true, false, false);
	
	vec params_nqs = {
					  sigma, sigm_2, sigm_4,
					  omg, omg_2,
					 };

    nqs n;
	n.N = N;
	n.M = M; 
	n.set_params(params_nqs);
    n.initialize();
	
	vec a_stored = colvec(M);
	vec b_stored = colvec(N);
	mat W_stored = mat(M, N);
	
	a_stored = n.a;
	b_stored = n.b;
	W_stored = n.W;
    
    for (double gamma: gamma_vals){
			n.a = a_stored;
			n.b = b_stored;
			n.W = W_stored;
      
            int n_sims = 50;
            int i = 0; 

            string base_filename = "weights.csv";
            
            while(i < n_sims){

                    retval result;

                    // Add 1000 to int to get proper sorting of filenames.
                    string filename = "/b_data/iteration_" 
                            + to_string(1000 + i)
                            + "_gamma_"
                            + to_string(gamma) 
                            + ".csv";
                    result = D.solve(n, filename);
                    
                    colvec a_update = colvec(M);
                    colvec b_update = colvec(N);
                    mat w_update = mat(M, N);
               
                    a_update = 2*(result.exp_vals.prod_E_grad_a 
                            - result.exp_vals.grad_a * result.el_exp);

                    b_update =  2*(result.exp_vals.prod_E_grad_b 
                            - result.exp_vals.grad_b * result.el_exp);

                    w_update =  2*(result.exp_vals.prod_E_grad_W 
                            - result.exp_vals.grad_W * result.el_exp);
                    
                    n.a += - gamma * a_update;
                    n.b += - gamma * b_update;
                    n.W += - gamma * w_update;
                    /*
                    cout << "WEIGTS" << endl;
                    n.a.print();
                    n.b.print();
                    n.W.print();
                    cout << "###########" << endl;
                    cout << "Iteration " << i << endl;
					*/
                    cout << "E_l = "<< result.el_exp << " | gamma " << gamma << endl;
                    timefile << result.time_spent << endl;
                    i ++;
        }// END gamma loop
    timefile.close();
    }
}
