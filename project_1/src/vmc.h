using namespace arma;
using namespace std;

class vmc {
    public:
        int n_dim, N, n_carlos;
        double step;
    
    private: 
        bool importance_sampling = false; 
        bool gradient_desc = false; 
        bool anal_derivative = false; 
        bool interactive = false;
        mat R;
        mat R_p;
        vector<vector<double>> simul_values;
	double best_alpha;
	double best_beta;
	double best_mean_E;
	double index_best;

    public :
        vmc(int n_dim, int N, int n_carlos, double step);
        void set_simul_params(bool imp, bool grad, bool anal, bool inter);
        void print_params();
        int solve(double alpha_start, double alpha_stop, double alpha_increment,
        double beta_start, double beta_stop, double beta_increment);
        void output_result();

    private :
        vector<double> carlo(double alpha, double beta);
        int variational_mc_naive(double beta_start, double beta_stop, double beta_increment,
                double alpha_start, double alpha_stop, double alpha_increment,
		 
                );

        void (*prop_move)();
        double local_energy_noninter(double alpha, double beta);
        double prob_dens_ratio_noninter(mat R_curr, double alpha, double beta);
        double prob_dens_ratio_inter(mat R_curr, double alpha, double beta);
        vector<double> output(mat E_l);
 


};
