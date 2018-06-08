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
    metadata(){}
    metadata(int M, int N){
        grad_a = colvec(M, fill::zeros);
        grad_b = colvec(N, fill::zeros);
        grad_W = mat(M, N, fill::zeros);

        prod_E_grad_a = colvec(M, fill::zeros);
        prod_E_grad_b = colvec(N, fill::zeros);
        prod_E_grad_W = mat(M, N, fill::zeros);
    }
    double *exp_E;
    colvec prod_E_grad_a;
    colvec grad_a;

    colvec grad_b;
    colvec prod_E_grad_b;

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
        bool compute_gibbs;
        int N_mc, N, M;
        int N_p, N_d;
        double sigma, omega, gamma;
   protected:
        colvec R;
        random_device rd;  //Will be used to obtain a seed for the random number engine
        mt19937 *gen; //Standard mersenne_twister_engine seeded with rd()
        colvec gradient_a, gradient_b;
		mat gradient_w;

    public:
        void monte_carlo(nqs& psi_t, metadata& exp_vals);
        retval solve(nqs& psi, string filename);
        void generate_positions(double step_int);
        void gradient_descent(nqs& psi_t, metadata& exp_vals, double E_l);
        void set_params(int N_p_in, int N_d_in, int N, int M, int mc_cycles,
                        bool obd_bool, bool meta_bool, bool gibbs_bool);
   protected:
        virtual double metropolis_hastings(nqs& psi_t, double prev_E_l) = 0;
};

class Importance: public vmc{
    protected:
        double metropolis_hastings(nqs& psi_t, double prev_E_l);
};

class NaiveMh: public vmc{

    protected:
        double metropolis_hastings(nqs& psi_t, double prev_E_l);
};

class Gibbs: public vmc{

    protected:
        double metropolis_hastings(nqs& psi_t, double prev_E_l);
};

