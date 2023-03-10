#include <stdio.h>
#include "WaveSolver.h"

int main(){

    double space_length = 2*M_PI;
    int N = 100;
    int Nt = 1;
    double space_step = space_length/(double(N)+1.);
    double time_step = space_step/2.;

    double** data = new double*[Nt];
    for (int i = 0; i < Nt; i++) data[i] = new double[N];

    WaveSolver wv;

    //giving initial data
    for (int j = 0; j < N; j++) data[0][j] = cos(j*space_step);

    for (int i = 0; i < N; i++){
        wv.TimeSimpleDiff1(data, Nt, N, space_step, time_step, wv.PFirstDerSpaceCenteredDiff2(data, Nt, N, space_step));
        Nt++;
    }
    
    wv.Write("graphics/output.txt", data, Nt-1, N, space_step, time_step);

    return 0;
}