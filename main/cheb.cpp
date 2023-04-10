#include <stdio.h>
#include "WaveSolver.h"
#include "Interpolator.h"

int main(){

    int N(15), Nt(1), res_inc(4);
    double space_length(1.), tf(5.);
    double space_step = space_length/(double(N));
    double time_step  = 0.001*space_step;
    int sizet = int(tf/time_step);

    double** data      = new double*[Nt];
    double** dot_data  = new double*[Nt];
    double*  chebx     = new double [N];
    double** precision = new double*[sizet+1];
    double*  prec_cheb = new double [res_inc*N];

    data[0]     = new double[N];
    dot_data[0] = new double[N];

    WaveSolver wv;
    wv.SetMethod("PSV"); //Sets the Pseudospectral as the chosen method

    for(int i = 0; i < sizet+1; i++)   precision[i] = new double[res_inc*N];
    for(int i = 0; i < res_inc*N; i++) prec_cheb[i] = cos(i*M_PI/(N*res_inc-1));

    //Initial data
    for (int j = 0; j < N; j++){
        chebx[j]       = cos(j*M_PI/(N-1));
        data[0][j]     = initialf(chebx[j]);
        dot_data[0][j] = initialdf(chebx[j]);
    }

    //Solving using PseudoSpectral method
    wv.TimeWaveRK4(data, dot_data, sizet, N, space_step, time_step);

    //Interpolating between the points given from wv
    Interpolator newton(chebx, data[0], N);

    for(int i = 0; i < sizet+1; i++){
        newton.SetY(data[i]);
        for(int j = 0; j < res_inc*N; j++) precision[i][j] = newton.Newton(cos(j*M_PI/(res_inc*N)));
    }

    //Writing the Output
    wv.Write("graphics/output.txt", precision, sizet+1, res_inc*N, time_step, prec_cheb);

    //Erasing everything
    for (int i = 0; i < sizet+1; i++){
        delete[] data[i]; delete[] dot_data[i]; delete[] precision[i];
    }

    delete[] data;  delete[] dot_data; delete[] precision;
    delete[] chebx; delete[] prec_cheb;

    return 0;
}