#include "../include/gaussian_noninter.h"
#include "../include/vmc.h"

int main(int argc, char *argv[]){
    GaussianNonInterNumeric g;
    
    importance D;
    D.set_params(1, 1, 1, 10);
    D.solve(&g, false);

    return 0;
}
