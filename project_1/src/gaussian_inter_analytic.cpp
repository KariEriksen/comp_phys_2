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

double GaussianInterAnalytic::E_l(mat R){
    double _psi = evaluate(R);
    double _laplace_psi, V_int = laplace(R);
    double V_ext = 0.5 * (double) as_scalar(accu(sum(square(R))));

    return - 0.5 * _laplace_psi +  V_ext + V_int;
}

double GaussianInterAnalytic::evaluate(mat R){
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

    vec rk, rj, ri;
    vec rkj, rki;

    double r_kj;
    double r_ki;

    double V_int = 0;

    for(int k = 0; k < N_p; k++){
        //Number of dimensions
        for(int dk = 0; dk < N_d; dk++){
            rk(k) = R(k, dk);
        }

        for(int j =  0; j < N_p; j++){

            //Number of dimensions
            for(int dj = 0; dj < N_d; dj++){
                rj(j) = R(j, dj);
            }

            rkj = rk - rj;

            r_kj = sqrt(sum(rkj*rkj));

            //r_kj = sqrt(xkj_sq + ykj_sq + zkj_sq);
            //r_kj = sqrt(sum(squared(R.row(j) - R.col(k)));
            //r_kj = D(k,j);

            if (j != k){

                sum_1 += ((sum(rkj))/r_kj)*(-a/(a*r_kj - r_kj*r_kj));

                part_1 = (a*(a-r_kj))/(r_kj*r_kj*(a - r_kj)*(a - r_kj));
                part_2 = 2.0/r_kj;
                part_3 = -a/(a*r_kj - r_kj*r_kj);

                sum_3 += part_1 + part_2 + part_3;
            }

            for(int i = 0; i < N_p; i++){

                //Number of dimensions
                for(int di = 0; di < N_d; di++){
                    ri(i) = R(i, di);
                }

                rki = rk - ri;

                r_ki = sqrt(sum(rki*rki));
                //r_ki = sqrt(xki_sq + yki_sq + zki_sq);
                //r_ki = D(k,i);

                if (i != k){
                    sum_2 += (-a/(a*r_ki - r_ki*r_ki))*(-a/(a*r_kj - r_kj*r_kj));
                }
            }

            if (k < j){
                if (r_kj < a){
                    V_int += 1;
                }
                else{
                    V_int += 0;
                }
            }
        }
    }

    double term1 = -4*alpha - 2*alpha*beta + 4*alpha_sq*(sum(rk*rk));
    double term2 = -4*alpha*(sum(rk))*sum_1;

    double scnd_der = term1 + term2 + sum_2 + sum_3;
    return scnd_der, V_int;
}

double GaussianInterAnalytic::drift_force(mat R){
    //Add this for the drift force to be used in importance sampling

    double sum_1 = 0;

    double alpha = params[0];
    double a = params[3];

    double xk, yk, zk;
    double xj, yj, zj;
    double xkj, ykj, zkj;
    double xkj_sq, ykj_sq, zkj_sq;

    double r_kj;
    double term = 0;

    for(int k = 0; k < N_p; k++){
        for(int j =  0; j < N_p; j++){

            xk = R(k,0);
            yk = R(k,1);
            zk = R(k,2);

            xj = R(j,0);
            yj = R(j,1);
            zj = R(j,2);

            xkj = (xk - xj);
            ykj = (yk - yj);
            zkj = (zk - zj);

            term = (-2*alpha*(xk - yk - zk));
            //r_kj = D(k,j);
            r_kj = sqrt(xkj_sq + ykj_sq + zkj_sq);

            if (j != k){
                sum_1 += (-a*(xkj + ykj + zkj))/(a*r_kj*r_kj - r_kj*r_kj*r_kj);

            }
        }
    }

    double first_der = term + sum_1;
    return first_der;
}


double GaussianInterAnalytic::ratio(mat R, mat R_p, int k){
    double eval_R = evaluate(R);
    double eval_R_p = evaluate(R_p);


    double prop = eval_R_p / eval_R;
    return prop*prop;
}

void GaussianInterAnalytic::update(){
    D = D_p;
}

void GaussianInterAnalytic::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    params = params_i;

}






