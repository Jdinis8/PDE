#include <stdio.h>
#include "WaveSolver.h"

int main(){
    double space_length = 1, tf = 2;
    int N(10);
    int Nt(1);
    int f(2), order(2);
    double space_step = space_length / (double(N));
    double time_step = 0.001 * space_step;
    int sizet = int(tf / time_step);

    double **data = new double *[Nt];
    double **dot_data = new double *[Nt];
    double *cheb = new double[N];

    for (int i = 0; i < Nt; i++)
    {
        data[i] = new double[N];
        dot_data[i] = new double[N];
    }

    for (int j = 0; j < N; j++)
    {
        data[0][j] = exp(-(cos((j * M_PI) / N)) * (cos((j * M_PI) / N)) / 0.1);
        dot_data[0][j] = 0.;
        cheb[j] = cos((j * M_PI) / N);
    }

    WaveSolver wv;

    std::vector<double> FD = wv.PSecondDerSpaceCenteredDiff2(data, 1, N, space_step);
    std::vector<double> PS = wv.PseudoSpectral(data, 1, N, cheb);
    std::cout << PS[0] << std::endl;

    for(int i = 0; i < int(FD.size()); i++){
        std::cout << FD[i] << " " << PS[i] << " " << FD[i]-PS[i] << std::endl;
    }

    for (int i = 0; i < sizet; i++)
    {
        delete[] data[i];
        delete[] dot_data[i];
    }

    delete[] data;
    delete[] dot_data;
    delete[] cheb;
    return 0;
}