#include <stdio.h>
#include "WaveSolver.h"

int main()
{

    double space_length = 2 * M_PI, x0 = -M_PI;
    int N = 100;
    int Nt = 1;
    int sizet = 10000;
    double space_step = space_length / (double(N));
    double time_step = 0.001;

    double **data = new double *[Nt];

    //here we are solving u' = dot(u)

    for (int i = 0; i < Nt; i++) data[i] = new double[N];
    WaveSolver wv;

    // giving initial data
    for (int j = 0; j < N; j++) data[0][j] = exp(-(x0 + j * space_step) * (x0 + j * space_step));

    for (int i = 0; i < sizet; i++)
    {
        wv.TimeSimpleDiff1(data, Nt, N, space_step, time_step, wv.PFirstDerSpaceCenteredDiff2(data, Nt, N, space_step));
        Nt++;
    }

    wv.Write("graphics/output.txt", data, Nt-1, N, space_step, time_step, x0);
    
    return 0;
}