#include "../include/nqs.h"
#include "../include/wavefunc.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

    int N, M;
    int N_d, N_mc, mc_exp;
    double gamma, omg, sigma, step;

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

    gamma = 0.01; omg = 1; sigma = 1; step = 0.1;
    double omg_2 = omg*omg;
    double sigm_2 = sigma*sigma;
    double sigm_4 = sigm_2*sigm_2;

    vec params_nqs = {
                      (double) M, (double) N,
                      sigma, sigm_2, sigm_4,
                      omg, omg_2,
                      gamma
                     };

    nqs n;
    NaiveMh D;

    D.step = step;
    n.set_params(params_nqs);

    vector<double> result;
    string filename = "filename";
    result = D.solve(&n, filename);


}
