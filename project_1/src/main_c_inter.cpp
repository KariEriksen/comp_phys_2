#include "../include/gaussian_inter_analytic.h"
#include "../include/wavefunc.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){
	// Derived object of WaveFunc class 

	double beta, step, h, a, dt; 
	int N_p, N_d, N_mc, mc_exp;
	beta = 1; step = 1; h = 1e-4; a = 0.0043;
	if( argc < 3){
		cout << "Wrong usage" << endl;
		cout << "Needs args: N_p N_d mc_exp dt" << endl;
		exit(1);
	}
	else{
		N_p = atoi(argv[1]); N_d = atoi(argv[2]); mc_exp = atoi(argv[3]); dt = atof(argv[4]);
	}
	N_mc = pow(2, mc_exp);

	Importance D;
	GaussianInterAnalytic g;

	vector<vector<double>> analytic_results; 

	double alpha_start, alpha_end, alpha_step;
	alpha_start = 0.1; alpha_end = 1; alpha_step = 0.05;
	int num_sims = (alpha_end- alpha_start)/alpha_step;
	double *alpha_array = new double[num_sims];

	for(int i = 0; i < num_sims ; i++){
		alpha_array[i] = alpha_start + i*alpha_step ;
	}

	for(int i = 0; i < num_sims; i++){
		double alpha = alpha_array[i];
		string sim_type_a = "IM_INA";
		vector<double> params = {alpha, alpha*alpha, beta, a, dt};
		g.set_params(params, N_d, N_p);
		//must be called or else you literally have no random numbers
		//args are alpha, beta, N_particles, N_dims, N_mccycles
		D.set_params(alpha, beta, N_p, N_d, N_mc);
		D.step = step;
		//reference must be passed or else you get static linking 
		//
		vector<double> result; 
		string filename = "../data/"+ sim_type_a+ "_a_" + to_string(alpha) + 
			"_b_" + to_string(beta) +
			"_step_" + to_string(step)+
			"_np_" + to_string(N_p)+
			"_nd_" + to_string(N_d)+
			"_dt_" + to_string(dt)+
			".csv";
		result = D.solve(&g, filename);
		analytic_results.push_back(result);
	}

	string meta_filename = "../data/IM_INA_meta_np_" + to_string(N_p)+
		"_nd_" + to_string(N_d)+   
		+"_data"
		+".csv";
	ofstream meta_file(meta_filename);

	meta_file << "alpha,analytic_energy,analytic_time,dt" << endl;
	double as = alpha_start;

	for(int i = 0; i<num_sims ; i++){
		meta_file << as << "," << analytic_results[i][0] << "," << analytic_results[i][1] << "," << dt << endl;
		//meta_file << "," << numeric_results[i][0] << "," << numeric_results[i][1] << endl;
		as += alpha_step;
	}

	meta_file.close();
	return 0;
}
