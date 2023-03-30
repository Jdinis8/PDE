#include <stdio.h>
#include "WaveSolver.h"

int main(){
    double space_length = 1, tf = 2;
    int N(21);
    int Nt(1);
    double space_step = space_length / (double(N));
    double time_step = 0.001 * space_step;
    int sizet = int(tf / time_step);

    double **data = new double *[Nt];
    double **dot_data = new double *[Nt];

    for (int i = 0; i < Nt; i++)
    {
        data[i] = new double[N];
        dot_data[i] = new double[N];
    }

    double x(0.);

    for (int j = 0; j < N; j++)
    {
        x = cos(j * M_PI / (N-1));
        data[0][j] = cos(x);
        dot_data[0][j] = 0.;
    }

    WaveSolver wv;

    std::vector<double> PS = wv.PseudoSpectral(data, 1, N);

    for(int j = 0; j < N; j++){
        x = cos(j*M_PI/(N-1));
        std::cout << PS[j] << " " << -cos(x) << std::endl;
    }


    for (int i = 0; i < 1; i++)
    {
        delete[] data[i];
        delete[] dot_data[i];
    }
    delete[] data;
    delete[] dot_data;
    return 0;
}