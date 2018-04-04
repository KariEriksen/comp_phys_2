#define arma_incl
#ifdef arma_incl
    #include <armadillo>
#endif

#include <random>
#include "../include/wavefunc.h"

using namespace std;
using namespace arma;

struct metadata{
    double *exp_E;
    double *prod_R; 
    double *prod_R_exp_E;
};


class vmc{
    public:
        double step;
        bool compute_extra; 
        double a, b;
        int N_p, N_mc, N_d;
   protected:
        mat R;
        random_device rd;  //Will be used to obtain a seed for the random number engine
        mt19937 *gen; //Standard mersenne_twister_engine seeded with rd()

    public:
        void monte_carlo(WaveFunc *psi_t, metadata *exp_vals);
        vector<double> solve(WaveFunc *psi, string filename);
        void generate_positions(double step_int);
        void set_params(double a, double b, int N, int dim,int mc_cycles, bool meta_bool = false);
    protected:
        virtual double metropolis_hastings(WaveFunc *psi_t, double prev_E_l) = 0;
};

class Importance: public vmc{
    protected:
        double metropolis_hastings(WaveFunc *psi_t, double prev_E_l);
};

class NaiveMh: public vmc{

    protected:
        double metropolis_hastings(WaveFunc *psi_t, double prev_E_l);
};

