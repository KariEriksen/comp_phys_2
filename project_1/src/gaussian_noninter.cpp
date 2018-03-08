#include "../include/gaussian_noninter.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianNonInterNumeric::GaussianNonInterNumeric() : WaveFunc(){}

mat GaussianNonInterNumeric::evaluate(mat R){
    double alpha = params[0];
    double beta = params[1];

    cout << "gauss_nonint evaluate " << endl;
    return exp(- alpha * square(R));
    
}
mat GaussianNonInterNumeric::laplace(mat R){
    return mat();
}

mat GaussianNonInterNumeric::nabla(mat R){
    return mat();
}
        

double GaussianNonInterNumeric::proportion(mat R, mat R_p){
    return 0;
}

void GaussianNonInterNumeric::set_params(vector<double> params, int N_d, int N_p){
    this -> N_d = N_d;
    this -> N_p = N_p;
    this -> params = params;
}
