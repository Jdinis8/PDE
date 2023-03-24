#include <stdio.h>
#include "WaveSolver.h"

int main(){

    double space_length = 1, x0 = 0, tf = 2;
    int N(200);
    int Nt(1);
    int f(2), order(2);
    double space_step = space_length / (double(N));
    double time_step = 0.001 * space_step;
    int sizet = int(tf / time_step);


    double **data = new double *[Nt];
    double **dot_data = new double *[Nt];
    double **copydata = new double *[Nt];
    double **copydot_data = new double *[Nt];

    for (int i = 0; i < Nt; i++){
        data[i] = new double[N];
        dot_data[i] = new double[N];
        copydata[i] = new double[N];
        copydot_data[i] = new double[N];
    }
}