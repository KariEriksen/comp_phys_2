#include "../include/gaussian_inter_analytic.h"
#include "../include/gaussian_inter_numeric.h"
#include "../include/gaussian_noninter_analytic.h"
#include "../include/gaussian_noninter_numeric.h"
#include "../include/vmc.h"

using namespace std;
using namespace arma;

//g++ -c gaussian_inter_analytic.cpp naive_mh.cpp vmc.cpp wavefunc.cpp -fPIC --std=c++11 -larmadillo
//g++ gaussian_inter_analytic.o naive_mh.o vmc.o wavefunc.o main_test_2.cpp -o main2 -fPIC -std=c++11 -larmadillo

int main(int argc, char *argv[]){
    // Derived object of WaveFunc class
    GaussianInterAnalytic g;
    //GaussianInterNumeric g;
    //GaussianNonInterAnalytic g;
    //GaussianNonInterNumeric g;

    NaiveMh D;
    //Importance D;
    double alpha, alpha_sq, beta, step, h, a;
    int N_p, N_d, N_mc;
    alpha = 0.5; alpha_sq = alpha*alpha;
    step = 0.1; h = 1e-5;
    N_p =5; N_d = 3; N_mc = 5;
    //beta = 2.82843; a = 0.00433;
    beta = 1.0; a = 0.0;

    //vector<double> params = {alpha, beta, h};
    //vector<double> params = {alpha, alpha_sq, beta};
    vector<double> params = {alpha, alpha_sq, beta, a};
    //vector<double> params = {alpha, beta, h, a};
    g.set_params(params, N_d, N_p);
    //must be called or else you literally have no random numbers
    //args are alpha, beta, N_particles, N_dims, N_mccycles
    D.set_params(alpha, beta, N_p, N_d, N_mc);
    D.step = step;
    //reference must be passed or else you get static linking
    //
    vector<double> result;
    string filename = "filename";
    result = D.solve(&g, filename);
    cout << "------------------------" << endl;
    cout << "E/N " << result[0]/N_p << endl;
    cout << "alpha " << alpha << endl;
    cout << "------------------------" << endl;

    return 0;
}
