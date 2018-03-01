#include <iostream>
#include <armadillo>
#include "vmc.h"

using namespace std;

int main(){
    /* Example usage 
    *vmc vmc_solver(1, 1, 1e2, 0.1);   
    *vmc_solver::solve("importance", "analytic", "none", false);
    *vmc_solver::output();
    */

    vmc solver(1, 1, 1e2, 0.1);
    solver.print_params();
    solver.solve(-1, 2, 0.001, 0, 1, 0.1);
    return 0;
}
