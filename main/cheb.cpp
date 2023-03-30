#include <stdio.h>
#include "WaveSolver.h"
#include "Interpolator.h"

int main(){

    double space_length = 1, x0 = 0, tf = 2;
    int N(15);
    int Nt(1);
    double space_step = space_length/(double(N));
    double time_step = 0.001*space_step;
    int sizet = int(tf/time_step);

    double** data = new double*[Nt];
    double** dot_data = new double*[Nt];
    double* chebx = new double[N];

    for (int i = 0; i < Nt; i++){
        chebx[i] = cos(i*M_PI/N);
        data[i] = new double[N];
        dot_data[i] = new double[N];
    }

    WaveSolver wv;
    wv.SetMethod("PSV");
    double x = 0;

    //giving initial data
    for (int j = 0; j < N; j++){
        x = cos((j * M_PI) / (N-1));
        data[0][j] = exp(-x * x / 0.1);
        dot_data[0][j] = 0.;
    }

    wv.TimeWaveRK4(data, dot_data, sizet, N, space_step, time_step);
 
    double** precision = new double*[sizet+1];

    for(int i = 0; i < sizet+1; i++) precision[i] = new double[4*N];

    Interpolator newton(chebx, data[0], N);
    for(int i = 0; i < sizet+1; i++){
        newton.SetY(data[i]);
        for(int j = 0; j < 4*N; j++) precision[i][j] = newton.Newton(cos(j*M_PI/(4.*N)));
    }
    wv.Write("graphics/output.txt", data, sizet+1, N, time_step, chebx);

    std::cout << "wow" << std::endl;

    for (int i = 0; i < sizet+1; i++){
        delete[] data[i];
        delete[] dot_data[i];
        delete[] precision[i];
    }

    delete[] data;
    delete[] dot_data;
    delete[] precision;
    delete[] chebx;
    return 0;
}