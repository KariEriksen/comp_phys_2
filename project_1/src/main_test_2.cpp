#include "../include/gaussian_noninter_numeric.h"
#include "../include/gaussian_noninter_analytic.h"
#include "../include/gaussian_inter_analytic.h"
#include "../include/vmc.h"
#include "../include/wavefunc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){
    // Derived object of WaveFunc class
    //GaussianNonInterAnalytic g;
    GaussianInterAnalytic g;

    NaiveMh D;
    double alpha, alpha_sq, beta, step, h, a;
    int N_p, N_d, N_mc;
    alpha = 0.5; alpha_sq = alpha*alpha;
    step = 0.1; h = 1e-5;
    N_p = 5; N_d = 3; N_mc = 5;
    //a = 0.0043; beta = 2.82843;
    a = 0.0; beta = 1.0;

    //vector<double> params = {alpha, alpha_sq, h};
    vector<double> params = {alpha, alpha_sq, beta, a};
    g.set_params(params, N_d, N_p);
    //must be called or else you literally have no random numbers
    //args are alpha, beta, N_particles, N_dims, N_mccycles
    D.set_params(alpha, beta, N_p, N_d, N_mc);
    D.step = step;
    //reference must be passed or else you get static linking
    //

    string filename = "test.csv";
    vector<double> result;
    result = D.solve(&g, filename);

    cout << "exp E/N" << "    " << result[0]/N_p << endl;

    return 0;
}
