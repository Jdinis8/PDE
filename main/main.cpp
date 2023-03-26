#include <stdio.h>
#include "WaveSolver.h"

int main(){

    double space_length = 2, x0 = -1, tf = 2;
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

    //giving initial data
    for (int j = 0; j < N; j++){
        data[0][j] = exp(-(x0 + j * space_step) * (x0 + j * space_step)/0.1);
        dot_data[0][j] = 0.;
        copydata[0][j] = data[0][j];
        copydot_data[0][j] = dot_data[0][j];
    }

    wv.TimeWaveRK4(data, dot_data, sizet, N, space_step, time_step);

    wv.Write("graphics/output.txt", data, sizet, N, space_step, time_step, x0);

    std::vector<std::vector<double>> conv = wv.ConvergenceTest(copydata, copydot_data, sizet, N, space_step, time_step, f, order);

    std::vector<double> l2 = wv.L2NormTime(copydata, copydot_data, sizet, N, space_step, time_step, f);

    wv.Write("graphics/conv.txt", conv, sizet, N/(f*f), space_step*f*f, time_step, x0);

    wv.Write("graphics/l2.txt", l2, time_step);

    for (int i = 0; i < sizet+1; i++){
        delete[] data[i];
        delete[] dot_data[i];
    }

    delete[] data;
    delete[] dot_data;
    return 0;
}