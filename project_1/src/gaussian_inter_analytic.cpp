#include "../include/gaussian_inter_analytic.h"
#include "../include/wavefunc.h"

#include <armadillo>

using namespace std;
using namespace arma;

GaussianInterAnalytic::GaussianInterAnalytic() : WaveFunc(){}

void GaussianInterAnalytic::initialize(mat R){
    double a = params[3];

    for(int i = 0; i < N_p; i ++ ){
        mat temp_outer = R.row(i);
        for(int j = (i+1) ; j < N_p; j++){
            mat temp_inner = R.row(j);
            mat temp = temp_inner - temp_outer;
            double dist = (double) sqrt(accu(square(temp)));
            D(i, j) = 1 - a/dist;
        }
    }
}
double GaussianInterAnalytic::eval_corr(mat R, int k = -1){
    double a = params[3];
    double ret_val = 1;

    if(k != -1){
        D_p = D;

        mat temp_outer = R.row(k);
        for(int i = 0;  k > i ; i++){
            mat temp_inner = R.row(i);
            mat temp = temp_inner - temp_outer;
            double dist = (double) sqrt(accu(square(temp)));
            if(dist > a){
                D_p(i, k) = 1 - a/dist;
            }
            else{
                D_p(i, k) = 0;
            }
        }

        for(int i = k; i < N_p ; i++){
            mat temp_inner = R.row(i);
            mat temp = temp_inner - temp_outer;
            double dist = (double) sqrt(accu(square(temp)));

            if(dist > a){
                D_p(k, i) = 1 - a/dist;
            }
            else{
                D_p(k, i) = 0;
            }
        }

        for(int i = 0; i < N_p; i++){
            for (int j = (i + 1); j < N_p ; j++){
                ret_val *= D_p(i, j);
            }
        }
        return ret_val;

    }
    else{
        for(int i = 0; i < N_p; i++){
            for (int j = (i + 1); j < N_p ; j++){
                ret_val *= D(i, j);
            }
        }
        return ret_val;
    }
}
double GaussianInterAnalytic::eval_g(mat R){
    double alpha = params[0];
    double beta = params[2];

    mat R_c(size(R));
    R_c = R;
    /*
	if(N_d > 2){
        R_c.col(2) *= beta;
    }*/
    double ret_val = 0;
    double internal = accu(sum(square(R_c)));

    ret_val = (double) as_scalar(exp(-alpha *(internal)));
    return ret_val;
}

double GaussianInterAnalytic::evaluate(mat R){
    double tmp = eval_corr(R);
    return eval_g(R)*tmp;
}

double GaussianInterAnalytic::E_l(mat R){
    double tmp = eval_corr(R);
    double _psi = eval_g(R) *tmp;
    double _laplace_psi = laplace(R);
    double V_ext = 0.5 * (double) as_scalar(accu(sum(square(R))));

    return - 0.5 * _laplace_psi +  V_ext;
}

double GaussianInterAnalytic::laplace(mat R){
    double alpha = params[0];
    double alpha_sq = params[1];
    double beta = params[2];
    double a = params[3];

    vec psi_d = zeros<vec>(N_d);
    vec rk = zeros<vec>(N_d);
    vec rj = zeros<vec>(N_d);
    vec ri = zeros<vec>(N_d);
    vec rkj = zeros<vec>(N_d);
    vec rki = zeros<vec>(N_d);

    double r_kj = 0;
    double r_ki = 0;
    double laplace_return = 0;

    for(int k = 0; k < N_p; k++){

        double sum_2 = 0;
        double sum_3 = 0;
        double psi_l = 0;
        double part_1 = 0;
        double part_2 = 0;
        double part_3 = 0;
        double part_4 = 0;
        double part_5 = 0;
        double part_6 = 0;

        vec sum_1 = zeros<vec>(N_d);
        //Number of dimensions
        for(int dk = 0; dk < N_d; dk++){
            rk(dk) = R(k, dk);
        }

        psi_l = -4*alpha - 2*alpha*beta + 4*alpha_sq*(sum(rk%rk));
        psi_d = -2*alpha*rk;

        for(int j = 0; j < N_p; j++){

            //Number of dimensions
            for(int dj = 0; dj < N_d; dj++){
                rj(dj) = R(j, dj);
            }

            if(j != k){

                rkj = rk - rj;
                r_kj = sqrt(sum(rkj%rkj));
                double d_u_rkj = -a/(a*r_kj - r_kj*r_kj);

                sum_1 += (rkj)/r_kj*(d_u_rkj);

                part_1 = (a*(a-2*r_kj))/(r_kj*r_kj*(a - r_kj)*(a - r_kj));
                part_2 = 2.0/r_kj;
                part_3 = -a/(a*r_kj - r_kj*r_kj);

                sum_3 += part_1 + part_2*part_3;



                for(int i = 0; i < N_p; i++){

                    //Number of dimensions
                    for(int di = 0; di < N_d; di++){
                        ri(di) = R(i, di);
                    }

                    if(i != k){

                        rki = rk - ri;
                        r_ki = sqrt(sum(rki%rki));
                        double d_u_rki = -a/(a*r_ki - r_ki*r_ki);
                        part_4 = d_u_rki;
                        part_5 = d_u_rkj;
                        part_6 = (sum(rki%rkj))/(r_ki*r_kj);

                        sum_2 += part_4*part_5*part_6;
                    }
                }
            }

        }
        laplace_return += psi_l + 2*sum(psi_d%sum_1) + sum_2 + sum_3;
    }

    return laplace_return;
}

mat GaussianInterAnalytic::drift_force(mat R){
    double alpha = params[0];
    double a = params[3];

    vec rk = zeros<vec>(N_d);
    vec rj = zeros<vec>(N_d);
    vec rkj = zeros<vec>(N_d);
    vec deri_phi_k = zeros<vec>(N_d);
    vec deri_u_k = zeros<vec>(N_d);

    double r_kj;

    for(int k = 0; k < N_p; k++){
        //Number of dimensions
        for(int dk = 0; dk < N_d; dk++){
            rk(dk) = R(k, dk);
        }

        for(int j =  0; j < N_p; j++){
            //Number of dimensions
            for(int dj = 0; dj < N_d; dj++){
                rk(dj) = R(j, dj);
            }

            rkj = rk - rj;
            r_kj = sqrt(sum(rkj%rkj));

            deri_phi_k = (-4*alpha*(rk));

            if (j != k){
                deri_u_k += -a/(a*r_kj*r_kj - r_kj*r_kj*r_kj)*rkj;

            }
        }
    }

    mat first_der = deri_phi_k + deri_u_k;
    return first_der;
}

double GaussianInterAnalytic::ratio(mat R, mat R_p, int k){
    double eval_R = eval_g(R)*eval_corr(R, k);
    double eval_R_p = eval_g(R_p)*eval_corr(R_p, k);

    double prob = (eval_R_p*eval_R_p)/(eval_R*eval_R);
    return prob;
}

void GaussianInterAnalytic::update(){
    //D_p.print();
    D = D_p;
}

void GaussianInterAnalytic::set_params(vector<double> params_i, int N_d_i, int N_p_i){
    N_d = N_d_i;
    N_p = N_p_i;
    params = params_i;

    D = mat(N_p, N_p, fill::zeros);
    D_p = mat(N_p, N_p, fill::zeros);
}






