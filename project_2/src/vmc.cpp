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

void vmc::monte_carlo(nqs& psi_t, metadata& exp_vals){

	if(compute_extra){
		exp_vals.exp_E[0] = psi_t.E_l(R);
	}
	if(compute_gibbs){
		exp_vals.exp_E[0] = psi_t.E_l_gibbs(R);
	}

	if(compute_extra || compute_gibbs){
		gradient_descent(psi_t, exp_vals, exp_vals.exp_E[0]);
	}

	// WRITING THINGS OUT BCUZ DAMN, GRAPH, U UGLY
	int acceptance = 0;
	//ofstream acceptance_file("../data/c_data/acceptance_file");

	for(int i = 1; i < N_mc; i++) {
		double tmp = metropolis_hastings(psi_t, exp_vals.exp_E[i-1]);
		exp_vals.exp_E[i] = tmp;
		if(exp_vals.exp_E[i] !=  exp_vals.exp_E[i-1]) acceptance++;

		if(compute_extra || compute_gibbs){
			gradient_descent(psi_t, exp_vals, tmp);
		}
	}
	cout << to_string(acceptance / (double) N_mc) << "\n";
	//acceptance_file << to_string(acceptance / (double) N_mc) << "\n";
	//acceptance_file.close();
}

void vmc::set_params(int N_p_in, int N_d_in,
		int N_in, int M_in, int mc_cycles,
		bool meta_bool,
		bool obd_bool, bool gibbs_bool){

	gradient_a = colvec(M_in);
	gradient_b = colvec(N_in);
	gradient_w = mat(M_in, N_in);

	N_mc = mc_cycles;
	R = colvec(M_in);

	M = M_in;
	N = N_in;
	N_p = N_p_in;
	N_d = N_d_in;

	gen = new mt19937(rd());
	compute_extra = meta_bool;
	compute_gibbs = gibbs_bool;
}

void vmc::generate_positions(double step_int){
	uniform_real_distribution<double> dis (-1, 1);

	for(int i = 0; i < M; i++){
		//gen must be dereferenced to use
		double tmp = dis(*gen) * step_int;
		R(i) = tmp;
	}
}


retval vmc::solve(nqs& psi_t, string filename){
	// -> is dereferencing and  member access to methods of class.
	generate_positions(step);

	/*
	 *Accepting a state in which a particle pair is closer than the permitted
	 * distance "a" is unphysical, and a waste of MC-cycles. We guarantee that the
	 * first state has a non-zero probability to occur and thus also guarantee that any
	 * proposed thate that has a < |r_i - r_j | for any particle pair is rejected
	 */


	metadata all_exp(M,N);
	all_exp.exp_E = new double[N_mc];

	int start_s = clock();
	monte_carlo(psi_t, all_exp);
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

	// Output to file
	ofstream output_file("../data/"+filename);
	arma_e_l.save(output_file, csv_ascii);
	output_file.close();

	retval r;
	r.exp_vals = all_exp;
	r.time_spent = time_spent;
	r.el_exp = accu(arma_e_l)/(double) N_mc ; 
	return r;
}

void vmc::gradient_descent(nqs& psi_t, metadata& exp_vals, double E_l){

	// Visible biases are vectors of length M.
	// Hidden nodes and corresponding hidden biases are vectors of length N
	// Weight matrix MxN
	// Number of visible nodes should correspond to the number of particles, N_p,
	// and the number of dimensions in the system, N_d, such that
	// M = N_p*N_d

	//int M = R.size(0);
	//int N = W.size(1);

	double sigma_squared = psi_t.sigma_2;

	// Gradient a
	gradient_a = R - psi_t.a;
	gradient_a /= sigma_squared;
	// Gradient b
	for (int k = 0; k < N; k++) {

		double sum_xi_wij = 0;
		for(int l = 0; l < M; l++){
			sum_xi_wij += R(l)*psi_t.W(l, k);
		}

		double inner_sum = sum_xi_wij/sigma_squared;
		gradient_b(k) = 1/(1 + exp(-psi_t.b(k) - inner_sum));
	}
	// Gradient w_kn
	// Concider changing summation indices in .tex for clarity
	for (int k = 0; k < M; k++) {
		for (int n = 0; n < N; n++) {

			double sum_xi_wij = 0;
			for(int l = 0; l < M; l++){
				sum_xi_wij += R(l)*psi_t.W(l, n);
			}

			double temp = sum_xi_wij/sigma_squared;
			gradient_w(k, n) = R(k)/(1 + exp(-psi_t.b(n) - temp));
		}
	}
	gradient_w /= sigma_squared;

	if(compute_extra){

		exp_vals.grad_a += gradient_a;
		exp_vals.grad_b += gradient_b;
		exp_vals.grad_W += gradient_w;

		exp_vals.prod_E_grad_a += E_l * gradient_a;
		exp_vals.prod_E_grad_b += E_l * gradient_b;
		exp_vals.prod_E_grad_W += E_l * gradient_w;
	}

	if(compute_gibbs){

		exp_vals.grad_a += 0.5*gradient_a;
		exp_vals.grad_b += 0.5*gradient_b;
		exp_vals.grad_W += 0.5*gradient_w;

		exp_vals.prod_E_grad_a += E_l * 0.5*gradient_a;
		exp_vals.prod_E_grad_b += E_l * 0.5*gradient_b;
		exp_vals.prod_E_grad_W += E_l * 0.5*gradient_w;
	}

}
