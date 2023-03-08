#include "WaveSolver.h"

WaveSolver::WaveSolver(){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
}

////////////////////////
////////////////////////
////////////////////////
//////TIME METHODS//////
////////////////////////
////////////////////////
////////////////////////

void WaveSolver::TimeSimpleDiff1(double** data, int size_t, int size_x, double space_step, double time_step, std::vector<double> RHS){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    
    double** new_data = new double* [size_t+1];
    double** dummy = data;


    for (int i = 0; i < size_t+1; i++) new_data[i] = new double[size_x];
    
    for (int i = 0; i < size_t; i++)
        for (int j = 0; j < size_x; j++)
            new_data[i][j] = data[i][j];

    for (int j = 0; j < size_x; j++)
        new_data[size_t+1][j] = new_data[size_t][j] + time_step*RHS[j];

    //Replacing data with the updated data, erasing the old array
    data = new_data;
    for(int i = 0; i < size_t; i++) delete[] dummy[i];
    delete[] dummy;
}

void WaveSolver::TimeCenteredDiff2(double **data, int size_t, int size_x, double space_step, double time_step, std::vector<double> RHS){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    double **new_data = new double *[size_t + 1];
    double **dummy = data;

    for (int i = 0; i < size_t + 1; i++)
        new_data[i] = new double[size_x];

    for (int i = 0; i < size_t; i++)
        for (int j = 0; j < size_x; j++)
            new_data[i][j] = data[i][j];

    for (int j = 0; j < size_x; j++)
        new_data[size_t + 1][j] = new_data[size_t-1][j] + 2*time_step * RHS[j];

    // Replacing data with the updated data, erasing the old array
    data = new_data;
    for (int i = 0; i < size_t; i++)
        delete[] dummy[i];
    delete[] dummy;
}

/////////////////////////
/////////////////////////
/////////////////////////
//////SPACE METHODS//////
/////////////////////////
/////////////////////////
/////////////////////////

std::vector<double> WaveSolver::PFirstDerSpaceCenteredDiff2(double **data, int size_t, int size_x, double space_step){
        
    //Below: size_t-1 means the most recent set of points, size_x-1 means the last point is equal to the point just before the first; even more below, when we are calculating the last derivative, the first term after the last one of our array must be equal to the initial one, because we have periodic solutions in this method
   
    std::vector<double> derivative;
    derivative.push_back((data[size_t - 1][size_x - 1] - data[size_t - 1][1]) / (2. * space_step));

    for (int i = 1; i < size_x-1; i++)
        derivative.push_back((data[size_t-1][i-1] - data[size_t - 1][i+1]) / (2. * space_step));

    derivative.push_back((data[size_t - 1][size_x - 1] - data[size_t - 1][0]) / (2. * space_step));

    return derivative;
}

std::vector<double> WaveSolver::PSecondDerSpaceCenteredDiff2(double **data, int size_t, int size_x, double space_step){

    std::vector<double> derivative;
    derivative.push_back((data[size_t - 1][size_x - 1] - data[size_t - 1][1]) / (2. * space_step));
}