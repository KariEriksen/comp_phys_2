#define arma_incl
#ifdef arma_incl
    #include <armadillo>
#endif

#include <random>
#include "../include/wavefunc.h"

using namespace std;
using namespace arma;

class vmc{
    public:
        double a, b;
        int N, mc_cycles;
   protected:
        mat R;
        random_device rd;  //Will be used to obtain a seed for the random number engine
        mt19937 *gen; //Standard mersenne_twister_engine seeded with rd()
        uniform_real_distribution<double> *dis; 

    public:
        vector<double> monte_carlo();
        void solve(wavefunc psi, bool import);
        void set_params(double a, double b, int N, int mc_cycles);
    protected:
        virtual mat metropolis_hastings(mat R) = 0;
};

class importance: public vmc{
    protected:
         mat metropolis_hastings(mat R);
};

class naive_mh: public vmc{
    protected:
        mat metropolis_hastings(mat R);
};

