#include "../include/nqs.h"
#include "../include/wavefunc.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

    int N, M;
    int N_d, N_mc, mc_exp;
    double gamma, omega, sigma, step;

    int N_p = 2;

    if( argc < 3){
        cout << "Wrong usage" << endl;
        exit(1);
    }
    else{
        N_d = atoi(argv[1]); mc_exp = atoi(argv[2]);
        N = atoi(argv[3]);
    }

    N_mc = pow(2, mc_exp);
    M = N_p*N_d;

    gamma = 0.01; omega = 1; sigma = 1; step = 0.1;

    nqs n;
    NaiveMh D;

    D.step = step;

    vector<double> result;
    string filename = "filename";
    result = D.solve(&n, filename);


}
