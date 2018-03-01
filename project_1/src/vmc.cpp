#include <armadillo>
#include <random>
#include "vmc.h"

using namespace std;
using namespace arma;

vmc::vmc(int n_dim, int N, int n_carlos, double step){
    n_dim = n_dim;
    N = N;
    n_carlos = n_carlos;
    step = step;
}
