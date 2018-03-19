#include "../include/gaussian_inter_numeric.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianInterNumeric::GaussianInterNumeric() : WaveFunc(){}

double GaussianInterNumeric::v_int(mat R){
    double a = params[3];

    for(int i = 0; i < N_p; i++){
        double tmp1 = (double) accu(sqrt(sum(square(R.row(i)))));
        for (int j = i + 1; j < N_p -1 ; j++){
            double tmp2 = (double) accu(sqrt(sum(square(R.row(j)))));
            if (std::abs(tmp1 - tmp2) < a) return 1e8;
        }
    }
    return 0;
}

double GaussianInterNumeric::jastrow(mat R){
    double a = params[3];
    double ret_val = 1;

    for(int i = 0; i < N_p; i++){
        double tmp1 = (double) accu(sqrt(sum(square(R.row(i)))));
        for (int j = i + 1; j < N_p -1 ; j++){
            double tmp2 = (double) accu(sqrt(sum(square(R.row(j)))));
            double comp = std::abs(tmp1 - tmp2);
            if (comp > a){
                ret_val *= 1 - a/comp;
            }
            else{
                return 0;
            }
        }
    }
    return ret_val;
}

double GaussianInterNumeric::E_l(mat R){
    double _psi = evaluate(R);
    double _laplace_psi = laplace(R);
    double V_ext = 0.5 * (double) as_scalar(accu(sum(square(R))));
    double V_int = v_int(R); 
    
    return - 0.5 * _laplace_psi/_psi +  V_ext + V_int;
}

double GaussianInterNumeric::evaluate(mat R){
    double alpha = params[0];
    double beta = params[1];

    mat R_c(size(R));
    R_c = R;
    if(N_d > 2){
        R_c.col(2) *= beta;
    }
    double ret_val = 0;
    double internal = accu(sum(square(R_c)));
    double jastrow_factors = jastrow(R);

    ret_val = (double) as_scalar(exp(-alpha *(internal)));
    return ret_val * jastrow_factors;
}
double GaussianInterNumeric::laplace(mat R){
    double h = params[2];
    double scnd_der = (evaluate(R-h) - 2* evaluate(R) + evaluate(R + h))/(h*h) ;
    return scnd_der;
}

double GaussianInterNumeric::nabla(mat R){
    double h = params[2];
    double der = (evaluate(R + h) - evaluate(R))/h;
    return der;
}


double GaussianInterNumeric::ratio(mat R, mat R_p){
    double eval_R = evaluate(R);
    double eval_R_p = evaluate(R_p);
    return eval_R / eval_R_p;
}

void GaussianInterNumeric::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    params = params_i;
}
