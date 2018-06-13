#include "../include/nqs.h"
#include "../include/wavefunc.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

    int N, M;
    int N_p, N_d, N_mc, mc_exp;
    double omg, sigma, step;

    N_p = 1;
    N_d = 1;
    mc_exp = 16;
    N = 3;

    N_mc = pow(2, mc_exp);
    M = N_p*N_d;
    
    colvec start_a = colvec(M);
    colvec start_b = colvec(N);
    mat start_w = mat(M, N);

    //double gamma_vals[] = {1e-3, 1e-2, 1e-1, 0.4};
    double gamma = 0.8;
    double sigma_vals[] = {0.9, 1, 1.3};
    ofstream timefile("../data/f_data/time_iter.csv");

    for (double sigma: sigma_vals){
            omg = 1; step = 0.45;
            double omg_2 = omg*omg;
            double sigm_2 = sigma*sigma;
            double sigm_4 = sigm_2*sigm_2;

            vec params_nqs = {
                              sigma, sigm_2, sigm_4,
                              omg, omg_2,
                              gamma
                             };
            nqs n;
            Gibbs D;

            D.step = step;
            D.set_params(N_p, N_d, N, M, N_mc, false, false, true);

            n.N = N;
            n.M = M;
            n.set_params(params_nqs);

            int n_sims = 100;
            int i = 0;

            if(i == 0){ 
                n.initialize();
                start_a = n.a;
                start_b = n.b;
                start_w = n.W;
            }
            else{
                n.a = start_a;
                n.b = start_b;
                n.W = start_w;
            }

            string base_filename = "weights.csv";

            while(i < n_sims){

                    retval result;

                    // Add 1000 to int to get proper sorting of filenames.
                    string filename = "../data/f_data/gibbs_iteration_"
                            + to_string(1000 + i)
                            + "_sigma_"
                            + to_string(sigma)
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
                    n.W.print();*/
                    cout << "###########" << endl;
                    cout << "Iteration " << i << endl;

                    cout << "E_l = "<< result.el_exp <<"\n" ;
                    //cout << endl;
                    cout << "sigma: " << sigma << endl;
                    timefile << result.time_spent << endl;
                    i ++;
        }// END gamma loop
    timefile.close();
    }
}
