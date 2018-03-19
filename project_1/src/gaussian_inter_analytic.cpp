#include "../include/gaussian_inter_analytic.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianNonInterAnalytic::GaussianNonInterAnalytic() : WaveFunc(){}

double GaussianInterAnalytic::v_int(mat R){
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

double GaussianInterAnalytic::E_l(mat R){
    double _psi = evaluate(R);
    double _laplace_psi = laplace(R);
    double V_ext = 0.5 * (double) as_scalar(accu(sum(square(R))));
    double V_int = v_int(R);

    return - 0.5 * _laplace_psi/_psi +  V_ext + V_int;
}

double GaussianInterAnalytic::evaluate(mat R){
    double alpha = params[0];
    double beta = params[1];

    mat R_c(size(R));
    R_c = R;
    if(N_d > 2){
        R_c.col(2) *= beta;
    }
    double ret_val = 0;
    double internal = accu(sum(square(R_c)));
    ret_val = (double) as_scalar(exp(-alpha *(internal)));
    return ret_val;
}

double GaussianInterAnalytic::laplace(mat R){

    double alpha = params[0];
    double alpha_sq = params[1];
    double beta = params[2];
    double a = params[3];

    double xk = R(k,j);
    double yk = R(k,j);
    double zk = R(k,j);

    double xk_sq = R(k,j);
    double yk_sq = R(k,j);
    double zk_sq = R(k,j);

    double sum_1 = 0;
    double sum_2 = 0;
    double sum_3 = 0;

    if (j != k){
        for(int k = 0; k < N_p; k++){
            for(int j =  0; j < N_d; j++){
                //NOT done!! and not correct
                r_kj = R(k,j);
                sum_1 += -a/(a*r_kj - r_kj*r_kj);

                part_1 = (a*(a-r_kj))/(r_kj*r_kj*(a - r_kj)*(a - r_kj));
                part_2 = 2.0/r_kj;
                part_3 = -a/(a*r_kj - r_kj*r_kj);
                sum_3 += part_1 + part_2 + part_3;
            }
        }
    }

    if (i,j != k){
        for(int k = 0; k < N_p; k++){
            for(int j =  0; j < N_d; j++){
                //NOT done!! and also not correct
                r_ki = R(k,i);
                r_kj = R(k,j);
                sum_2 += (-a/(a*r_ki - r_ki*r_ki))*(-a/(a*r_kj - r_kj*r_kj));
            }
        }
    }


    double term1 = -4*alpha - 2*alpha*beta + 4*alpha_sq*(xk_sq + yk_sq + zk_sq);
    double term2 = -4*alpha*(xk + yk + zk)*sum_1;

    double scnd_der = term1 + term2 + sum_2 + sum_3;
    return scnd_der;
}

double GaussianInterAnalytic::nabla(mat R){
    //Add this for the drift force to be used in importance sampling
    return 0;
}


double GaussianInterAnalytic::ratio(mat R, mat R_p){
    double eval_R = evaluate(R);
    double eval_R_p = evaluate(R_p);

    double prop = eval_R / eval_R_p;
    return prop;
}

void GaussianInterAnalytic::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    params = params_i;

}






