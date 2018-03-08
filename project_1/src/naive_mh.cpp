#include <armadillo>
#include "../include/vmc.h"

using namespace arma;

mat naive_mh::metropolis_hastings(mat R, WaveFunc *psi_t){
    return psi_t -> evaluate(R);
}

