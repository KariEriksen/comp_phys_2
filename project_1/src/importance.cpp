#include <armadillo>
#include "../include/vmc.h"

using namespace arma;

mat importance::metropolis_hastings(mat R, WaveFunc *psi_t){
    cout << "AGAIN, HELLO" << endl;
    psi_t -> evaluate(R);
    return R;
}
