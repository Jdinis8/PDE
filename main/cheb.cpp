#include <stdio.h>
#include "WaveSolver.h"

int main(){

    double space_length = 1, x0 = 0, tf = 2;
    int N(200);
    int Nt(1);
    int f(2), order(2);
    double space_step = space_length/(double(N));
    double time_step = 0.001*space_step;
    int sizet = int(tf/time_step);

    double** data = new double*[Nt];
    double** dot_data = new double*[Nt];
    double** copydata = new double*[Nt];
    double** copydot_data = new double*[Nt];

    for (int i = 0; i < Nt; i++){
        data[i] = new double[N];
        dot_data[i] = new double[N];
        copydata[i] = new double[N];
        copydot_data[i] = new double[N];
    }

    WaveSolver wv;
    wv.SetMethod("PS");

    //giving initial data
    for (int j = 0; j < N; j++){
        data[0][j] = exp(-(x0 + cos((j * M_PI) / N)) * (x0 + cos((j * M_PI) / N)) / 0.1);
        dot_data[0][j] = 0.;
        copydata[0][j] = data[0][j];
        copydot_data[0][j] = dot_data[0][j];
    }

    wv.TimeWaveRK4(data, dot_data, sizet, N, space_step, time_step);
    
    wv.Write("graphics/output.txt", data, sizet, N, space_step, time_step, x0);

    for (int i = 0; i < sizet+1; i++){
        delete[] data[i];
        delete[] dot_data[i];
    }

    delete[] data;
    delete[] dot_data;
    return 0;
}