#include "../include/gaussian_inter_analytic.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianInterAnalytic::GaussianInterAnalytic() : WaveFunc(){}

double GaussianInterAnalytic::E_l(mat R){
    double _over_psi = evaluate(R);
    double _laplace_psi = laplace(R);

    return _over_psi*_laplace_psi;

}
