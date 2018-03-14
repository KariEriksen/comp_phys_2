#define arma_incl
#ifdef arma_incl
    #include <armadillo>
#endif

#include <random>
#include "../include/wavefunc.h"
#include "../include/gaussian_noninter_numeric.h"
#include "../include/gaussian_noninter_analytic.h"
#include "../include/gaussian_inter_numeric.h"
#include "../include/gaussian_inter_analytic.h"

using namespace std;
using namespace arma;

class vmc{
    public:
        double step;       
        double a, b;
        int N_p, N_mc, N_d;
   protected:
        mat R;
        random_device rd;  //Will be used to obtain a seed for the random number engine
        mt19937 *gen; //Standard mersenne_twister_engine seeded with rd()

    public:
        vector<double> monte_carlo(WaveFunc *psi_t);
        vector<double> solve(WaveFunc *psi);
        void set_params(double a, double b, int N, int dim,int mc_cycles);
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

