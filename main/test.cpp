#include <stdio.h>
#include "WaveSolver.h"
#include <iomanip>
#include <limits>

int main(){
    double space_length = 1, tf = 2;
    int N(20);
    int Nt(1);
    double space_step = space_length / (double(N));
    double time_step = 0.001 * space_step;
    int sizet = int(tf / time_step);
    constexpr auto max_precision{std::numeric_limits<long double>::digits10 + 1};

    std::cout << std::setprecision(max_precision);

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
        data[0][j] = exp(x)*sin(5*x);
        dot_data[0][j] = 0.;
    }

    WaveSolver wv;

    //std::vector<double> PS = wv.PseudoSpectral(data, 1, N);
    std::vector<double> FFT = wv.FFTPseudoSpectral(data, 1, N);

    for(int j = 0; j < N; j++){
        x = cos(j*M_PI/(N-1));
        std::cout << cos(j*M_PI/(N-1)) << " " << FFT[j] << " " << exp(x) * (sin(5 * x) + 5 * cos(5 * x)) << " " << FFT[j] - exp(x) * (sin(5 * x) + 5 * cos(5 * x))  << std::endl;
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