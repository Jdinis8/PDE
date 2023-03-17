#include <stdio.h>
#include "WaveSolver.h"

int main()
{
    double space_length = 2 * M_PI, x0 = -M_PI, tf = 3;
    int N(200);
    int f(2), order(2);
    double space_step = space_length / (double(N));
    double time_step = 0.001 * space_step;
    int sizet = int(tf / time_step);

    WaveSolver wv;

    // we'll do several cases for calculating l2norm
    int iter(100);
    int N2(N);
    double sstep(space_step);
    double tstep(sstep);
    int sizet2(sizet);
    double **errors = new double *[2];
    errors[0] = new double[iter];
    errors[1] = new double[iter];

    for (int i = 0; i < iter; i++)
    {
        N2 = N + i * f * f;
        sstep = space_length / double(N2);
        tstep = 0.001 * sstep;
        sizet2 = int(tf / tstep);

        double **data2 = new double *[1];
        double **dot_data2 = new double *[1];
        data2[0] = new double[N2];
        dot_data2[0] = new double[N2];

        for (int j = 0; j < N2; j++)
        {
            data2[0][j] = exp(-(x0 + j * sstep) * (x0 + j * sstep));
            dot_data2[0][j] = 0.;
        }

        errors[0][i] = sstep;
        errors[1][i] = wv.L2NormStep(data2, dot_data2, sizet2, N2, sstep, tstep, f, order);
        std::cout << "h: " << sstep << " error l2norm: " << errors[1][i] << std::endl;

        delete[] data2[0];
        delete[] dot_data2[0];
        delete[] data2;
        delete[] dot_data2;
    }

    wv.Write("graphics/l2.txt", errors, iter);
}