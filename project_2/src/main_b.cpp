#include "../include/nqs.h"
#include "../include/wavefunc.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

    int N, M;
    int N_p, N_d, N_mc;
    double gamma, omega, sigma, step;

    //N = 3; N_p = 1; N_d = 2;
    N = atoi(argv[1]) ; N_p = atoi(argv[1]);
    N_d = atoi(argv[2]); mc_exp = atoi(argv[3]);

    M = N_p*N_d;

    N_mc = pow(2, mc_exp);

    gamma = 1; omega = 1; sigma = 1; step = 0.1;

    nqs n;
    NaiveMh D;

    D.step = step;

    vector<double> result;
    string filename = "filename";
    result = D.solve(&n, filename);


}
