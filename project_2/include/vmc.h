#define arma_incl
#ifdef arma_incl
    #include <armadillo>
#endif

#include <random>
#include "../include/wavefunc.h"
#include "../include/nqs.h"

using namespace std;
using namespace arma;

struct metadata{
    double *exp_E;
    mat prod_E_grad_a;
    mat grad_a;

    mat grad_b;
    mat prod_E_grad_b;

    mat grad_W;
    mat prod_E_grad_W;

    long int *obd;
};

struct retval{
    metadata exp_vals;
    double time_spent; 
    double el_exp;
};


class vmc{
    public:
		double dt;
        double step;
        bool compute_extra; 
        bool compute_obd;
        int N_mc, N, M;
        double sigma, omega, gamma;
   protected:
        int obd_n_bins;
        double bin_length;
        mat R;
        random_device rd;  //Will be used to obtain a seed for the random number engine
        mt19937 *gen; //Standard mersenne_twister_engine seeded with rd()
        double* obd_bins;
        void count_obd(mat radii, metadata* all_exp);
        mat gradient_a, gradient_b, gradient_w;

    public:
        void monte_carlo(nqs *psi_t, metadata *exp_vals);
        retval solve(nqs *psi, string filename);
        void generate_positions(double step_int);
        void gradient_descent(nqs *psi_t, metadata *exp_vals, double E_l);
        void set_params(int N, int M, int mc_cycles, bool obd_bool, bool meta_bool);
    protected:
        virtual double metropolis_hastings(nqs *psi_t, double prev_E_l) = 0;
};

class Importance: public vmc{
    protected:
        double metropolis_hastings(nqs *psi_t, double prev_E_l);
};

class NaiveMh: public vmc{

    protected:
        double metropolis_hastings(nqs *psi_t, double prev_E_l);
};

