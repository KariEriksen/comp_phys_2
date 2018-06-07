#include "../include/nqs.h"
#include "../include/wavefunc.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

    int N, M;
    int N_d, N_mc, mc_exp;
    double gamma, omg, sigma, step;

    int N_p = 1;

    /*
    if( argc < 3){
        cout << "Wrong usage" << endl;
        exit(1);
    }
    else{
        N_d = atoi(argv[1]); mc_exp = atoi(argv[2]);
        N = atoi(argv[3]);
    }
    */
    N_d = 1;
    mc_exp = 15;
    N = N_p * N_d;

    N_mc = pow(2, mc_exp);
    M = N_p*N_d;

    gamma = 1e-1; omg = 1; sigma = 0.9; step = 0.1;
    double omg_2 = omg*omg;
    double sigm_2 = sigma*sigma;
    double sigm_4 = sigm_2*sigm_2;
    
    vec params_nqs = {
                      sigma, sigm_2, sigm_4,
                      omg, omg_2,
                      gamma
                     };
    nqs n;
    NaiveMh D;

    D.step = step;
    D.set_params(N_p, N_d, N, M, N_mc, 1, 0);

    n.N = N;
    n.M = M; 
    n.set_params(params_nqs);
    n.initialize();
  
    int n_sims = 5; 
    int i = 0; 

    string base_filename = "weights.csv";

    while(i < n_sims){

        retval result;
        string filename = "/b_data/iteration_"+to_string(i)+".csv";
        result = D.solve(&n, filename);
        
        colvec a_update = colvec(M);
        colvec b_update = colvec(N);
        mat w_update = mat(M, N);
       
        a_update = 2*(result.exp_vals.prod_E_grad_a 
            - result.exp_vals.grad_a * result.el_exp);

        b_update =  2*(result.exp_vals.prod_E_grad_b 
            - result.exp_vals.grad_b * result.el_exp);

        w_update =  2*(result.exp_vals.prod_E_grad_W 
            - result.exp_vals.grad_W * result.el_exp);
 

        n.a = n.a - gamma * a_update;
        n.b = n.b - gamma * b_update;
        n.W = n.W - gamma * w_update;
        
        cout << "Iteration " << i << endl;
        cout << "E_l = "<< result.el_exp << endl;
        i ++;
    }
}
