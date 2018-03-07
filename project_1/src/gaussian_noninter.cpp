#include "../include/gaussian_noninter.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianNonInterNumeric::GaussianNonInterNumeric() : wavefunc(){}

mat GaussianNonInterNumeric::evaluate(mat R){
    double alpha = params[0];

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

void GaussianNonInterNumeric::set_params(vector<double> params){
    this -> params = params;
}
