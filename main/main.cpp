#include <stdio.h>
#include "WaveSolver.h"

int main(){

    double space_length = 2*M_PI;
    int N = 2000;
    int Nt = 1;
    int sizet = 10000;
    double space_step = space_length/(double(N));
    double time_step = 0.001;

    double** data = new double*[Nt];
    double** dot_data = new double*[Nt];

    for (int i = 0; i < Nt; i++){
        data[i] = new double[N];
        dot_data[i] = new double[N];
    }

    WaveSolver wv;

    //giving initial data
    for (int j = 0; j < N; j++){
        data[0][j] = cos(j*space_step);
        dot_data[0][j] = sin(2.*j * space_step);
    }

    for (int i = 0; i < N; i++){
        //wv.TimeSimpleDiff1(data, Nt, N, space_step, time_step, wv.PFirstDerSpaceCenteredDiff2(data, Nt, N, space_step));
        Nt++;
    }

    // wv.Write("graphics/output.txt", data, Nt-1, N, space_step, time_step);

    wv.TimeWaveRK4(data, dot_data, sizet, N, space_step, time_step);
    wv.Write("graphics/output.txt", data, sizet, N, space_step, time_step);

    for (int i = 0; i < sizet+1; i++)
    {
        delete[] data[i];
        delete[] dot_data[i];
    }

    delete[] data;
    delete[] dot_data;
    return 0;
}