#include "../include/gaussian_inter_analytic.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianInterAnalytic::GaussianInterAnalytic() : WaveFunc(){}

void GaussianInterAnalytic::initialize(mat R){
    double a = params[3];
    D = mat(N_p, N_p);
    D.zeros();
    for(int i = 0; i < N_p; i ++ ){
        mat temp_outer = R.row(i);
        for(int j = (i+1) ; j < N_p; j++){
            mat temp_inner = R.row(j);
            D(i, j) = 1 - a/(std::abs((double) accu(sum(temp_outer - temp_inner))));
        }
    }
}
double GaussianInterAnalytic::eval_corr(mat R, int k = -1){
    double a = params[3];
    double ret_val = 1;
    mat D_c = D;

    if(k != -1){
        mat r_k = R.row(k);
        for(int i = 0;  k > i ; i++){
            double r_ik = std::abs(sum(R.row(i) - r_k));
            if(r_ik > a){
                D(i, k) = r_ik;
            }
            else{
                D(i, k) = 0;
            }
        }
        for(int i = k; i < N_p ; i++){
            double r_ik = std::abs(sum(R.row(i) - r_k));
            if(r_ik > a){
                D(k, i) = r_ik;
            }
            else{
                D(k, i) = 0;
            }
        }
    }
    for(int i = 0; i < N_p; i++){
        for (int j = (i + 1); j < N_p ; j++){
            ret_val *= D(i, j);
        }
    }
    if(k != -1){
        D_p = D;
        D = D_c;
    }
    return ret_val;
}

double GaussianInterAnalytic::eval_g(mat R){
    double alpha = params[0];
    double beta = params[2];

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

double GaussianInterAnalytic::evaluate(mat R){
    return eval_g(R)*eval_corr(R);
}

double GaussianInterAnalytic::E_l(mat R){
    double _psi = evaluate(R);
    double _laplace_psi, V_int = laplace(R);
    double V_ext = 0.5 * (double) as_scalar(accu(sum(square(R))));

    return - 0.5 * _laplace_psi +  V_ext + V_int;
}

double GaussianInterAnalytic::laplace(mat R){

    double alpha = params[0];
    double alpha_sq = params[1];
    double beta = params[2];
    double a = params[3];

    double sum_1 = 0;
    double sum_2 = 0;
    double sum_3 = 0;

    double part_1 = 0;
    double part_2 = 0;
    double part_3 = 0;

    vec rk = zeros<vec>(N_d);
    vec rj = zeros<vec>(N_d);
    vec ri = zeros<vec>(N_d);
    vec rkj = zeros<vec>(N_d);
    vec rki = zeros<vec>(N_d);

    double r_kj;
    double r_ki;

    double V_int = 0;

    for(int k = 0; k < N_p; k++){

        //Number of dimensions
        for(int dk = 0; dk < N_d; dk++){
            rk(dk) = R(k, dk);
        }

        for(int j = 0; j < N_p; j++){

            //Number of dimensions
            for(int dj = 0; dj < N_d; dj++){
                rj(dj) = R(j, dj);
            }

            if(j != k){

                rkj = rk - rj;
                r_kj = sqrt(sum(rkj%rkj));

                sum_1 += -4*alpha*(sum(rk%rkj)/r_kj)*(-a/(a*r_kj - r_kj*r_kj));

                part_1 = (a*(a-r_kj))/(r_kj*r_kj*(a - r_kj)*(a - r_kj));
                part_2 = 2.0/r_kj;
                part_3 = -a/(a*r_kj - r_kj*r_kj);

                sum_3 += part_1 + part_2 + part_3;

            }

            for(int i = 0; i < N_p; i++){

                //Number of dimensions
                for(int di = 0; di < N_d; di++){
                    ri(di) = R(i, di);
                }

                if(i != k){

                    rki = rk - ri;
                    r_ki = sqrt(sum(rki%rki));

                    sum_2 += (sum(rki%rkj))/(r_ki*r_kj)*(-a/(a*r_ki - r_ki*r_ki))*(-a/(a*r_kj - r_kj*r_kj));
                }
            }

            if (k < j){
                if (r_kj <= a){
                    V_int += 1e10;
                }
                else{
                    V_int += 0;
                }
            }
        }
    }

    double term1 = -4*alpha - 2*alpha*beta + 4*alpha_sq*(sum(rk%rk));

    double scnd_der = term1 + sum_1 + sum_2 + sum_3;
    return scnd_der, V_int;
}

double GaussianInterAnalytic::drift_force(mat R){
    //Add this for the drift force to be used in importance sampling

    double sum_1 = 0;

    double alpha = params[0];
    double a = params[3];

    vec rk = zeros<vec>(N_d);
    vec rj = zeros<vec>(N_d);
    vec rkj = zeros<vec>(N_d);

    double r_kj;
    double term = 0;
    double addt = 0;

    for(int k = 0; k < N_p; k++){
        //Number of dimensions
        for(int dk = 0; dk < N_d; dk++){
            rk(dk) = R(k, dk);
            addt += rk(dk);
        }

        for(int j =  0; j < N_p; j++){
            //Number of dimensions
            for(int dj = 0; dj < N_d; dj++){
                rk(dj) = R(j, dj);
            }

            rkj = rk - rj;
            r_kj = sqrt(sum(rkj%rkj));

            term = (-4*alpha*(addt));
            //r_kj = D(k,j);

            if (j != k){
                sum_1 += (-a*(sum(rkj)))/(a*r_kj*r_kj - r_kj*r_kj*r_kj);

            }
        }
    }

    double first_der = term + sum_1;
    return first_der;
}

double GaussianInterAnalytic::ratio(mat R, mat R_p, int k){
    double eval_R = evaluate(R)*eval_corr(R, k);
    double eval_R_p = evaluate(R_p)*eval_corr(R_p, k);


    double prob = (eval_R_p*eval_R_p)/(eval_R*eval_R);
    return prob;
}

void GaussianInterAnalytic::update(){
    D = D_p;
}

void GaussianInterAnalytic::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    params = params_i;

}






