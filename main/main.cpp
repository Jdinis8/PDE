#include <stdio.h>
#include "WaveSolver.h"

int main(){

    double space_length = 2*M_PI, x0 = -M_PI;
    int N = 100;
    int Nt = 1;
    int sizet = 10000;
    double space_step = space_length/(double(N));
    double time_step = 0.001;

    double** data = new double*[Nt];
    double** dot_data = new double*[Nt];

    for (int i = 0; i < N; i++){
        data[i] = new double[N];
        dot_data[i] = new double[N];
    }

    WaveSolver wv;

    //giving initial data
    for (int j = 0; j < N; j++){
        data[0][j] = exp(-(x0 + j * space_step) * (x0 + j * space_step));
        dot_data[0][j] = 0.;
    }

    wv.TimeWaveRK4(data, dot_data, sizet, N, space_step, time_step);
    wv.Write("graphics/output.txt", data, sizet, N, space_step, time_step, x0);

    for (int i = 0; i < sizet+1; i++)
    {
        delete[] data[i];
        delete[] dot_data[i];
    }

    delete[] data;
    delete[] dot_data;
    return 0;
}