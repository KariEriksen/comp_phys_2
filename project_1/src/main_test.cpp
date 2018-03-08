#include "../include/gaussian_noninter.h"
#include "../include/vmc.h"

int main(int argc, char *argv[]){
    // Derived object of WaveFunc class 
    GaussianNonInterNumeric g;
    
    importance D;
    //must be called or else you literally have no random numbers
    D.set_params(1, 1, 1, 10);
    //reference must be passed or else you get static linking 
    D.solve(&g, false);

    return 0;
}
