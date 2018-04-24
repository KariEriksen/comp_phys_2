#include <armadillo>
#include "../include/vmc.h"
#include <random>

using namespace arma;
using namespace std;

double NaiveMh::metropolis_hastings(WaveFunc *psi_t, double prev_E_l){
    mat R_p(size(R));
    R_p = R;
	double beta = psi_t -> params[2]

    uniform_int_distribution<int> dis_r(0, N_p - 1);
    uniform_real_distribution<double> dis_step(-1, 1);
    uniform_real_distribution<double> dis_p(0, 1);

    int j = dis_r(*gen);
    for(int i = 0; i < N_d; i++){
        R_p(j, i) += dis_step(*gen) * step ;
    }
	if(N_d == 3){
		R_p(j,2) *= beta;
	}
    double eps = dis_p(*gen);
    double P = psi_t -> ratio(R, R_p, j);
    /*
    cout << "eps: " << eps << " | P: " << P << endl;
    cout << "R ------------------" << endl;
    R.print();
    cout << "R_p ----------------" << endl;
    R_p.print();
    */
    if(eps < P){
        R = R_p;
        psi_t -> update();
        return psi_t -> E_l(R);
    }
    else{
        return prev_E_l;
    }

}


