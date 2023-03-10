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

void WaveSolver::TimeSimpleDiff1(double**& data, int size_t, int size_x, double space_step, double time_step, std::vector<double> RHS){
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
        new_data[size_t][j] = new_data[size_t-1][j] + time_step*RHS[j];

    data = new_data;
    //Replacing data with the updated data, erasing the old array
    for(int i = 0; i < size_t; i++) delete[] dummy[i];
    delete[] dummy;
}

void WaveSolver::TimeCenteredDiff2(double**& data, int size_t, int size_x, double space_step, double time_step, std::vector<double> RHS){
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
        new_data[size_t][j] = new_data[size_t-1][j] + 2*time_step * RHS[j];

    // Replacing data with the updated data, erasing the old array
    data = new_data;

    for (int i = 0; i < size_t; i++)
        delete[] dummy[i];
    delete[] dummy;
}

void WaveSolver::TimeRK4(double**& data, int size_t, int size_x, double space_step, double time_step){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
}

/////////////////////////
/////////////////////////
/////////////////////////
//////SPACE METHODS//////
/////////////////////////
/////////////////////////
/////////////////////////
//////PERIODIC ONES//////
/////////////////////////

std::vector<double> WaveSolver::PFirstDerSpaceCenteredDiff2(double ** data, int size_t, int size_x, double space_step){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    //Below: size_t-1 means the most recent set of points, size_x-1 means the last point is equal to the point just before the first; even more below, when we are calculating the last derivative, the first term after the last one of our array must be equal to the initial one, because we have periodic solutions in this method
   
    std::vector<double> derivative;
    derivative.push_back((data[size_t - 1][size_x - 1] - data[size_t - 1][1]) / (2. * space_step));

    for (int i = 1; i < size_x-1; i++)
        derivative.push_back((data[size_t-1][i-1] - data[size_t - 1][i+1]) / (2. * space_step));

    derivative.push_back((data[size_t - 1][size_x - 2] - data[size_t - 1][0]) / (2. * space_step));

    return derivative;
}

std::vector<double> WaveSolver::PFirstDerSpaceCenteredDiff4(double **data, int size_t, int size_x, double space_step){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    //fourth order difference
    std::vector<double> derivative;

    if (size_x < 5){
        std::cout << "[WARNING] The size of your data is not enough to compute this approximation, please try another one!" << std::endl;
        return derivative;
    }

    //-(-2) + 8*(-1) - 8*(+1) + (+2)
    derivative.push_back((-data[size_t-1][2] + 8.*data[size_t-1][1] - 8.*data[size_t-1][size_x-1] + data[size_t-1][size_x-2]) / (12. * space_step));
    derivative.push_back((-data[size_t-1][3] + 8.*data[size_t-1][2] - 8.*data[size_t-1][0] + data[size_t-1][size_x-1]) / (12. * space_step));

    for (int i = 2; i < size_x - 3; i++)
        derivative.push_back((-data[size_t - 1][i+2] + 8. * data[size_t - 1][i+1] - 8. * data[size_t - 1][i-1] + data[size_t - 1][i-2]) / (12. * space_step));

    //for size_x-2 step
    derivative.push_back((-data[size_t - 1][0] + 8. * data[size_t - 1][size_x-1] - 8. * data[size_t - 1][size_x-3] + data[size_t - 1][size_x-4]) / (12. * space_step));
    derivative.push_back((-data[size_t - 1][1] + 8. * data[size_t - 1][0] - 8. * data[size_t - 1][size_x-2] + data[size_t - 1][size_x-3]) / (12. * space_step));

    return derivative;
}

std::vector<double> WaveSolver::PSecondDerSpaceCenteredDiff2(double **data, int size_t, int size_x, double space_step){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    //see explanation above for this first calculation and the last one

    std::vector<double> derivative;
    derivative.push_back((data[size_t - 1][size_x - 1] - 2.*data[size_t - 1][0] + data[size_t-1][1]) / (space_step * space_step));

    for (int i = 1; i < size_x-1; i++)
        derivative.push_back((data[size_t - 1][i-1] - 2. * data[size_t - 1][i] + data[size_t - 1][i+1]) / (space_step * space_step));

    derivative.push_back((data[size_t - 1][size_x - 2] - 2. * data[size_t - 1][size_x-1] + data[size_t - 1][0]) / (space_step * space_step));

    return derivative;
}

/////////////////////////
////NON-PERIODIC ONES////
/////////////////////////

std::vector<double> WaveSolver::FirstDerSpaceCenteredDiff2(double **data, int size_t, int size_x, double space_step){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    //we'll use forward difference for the first step and backward difference for the last one
    
    std::vector<double> derivative;

    //forward: -3*now + 4*after - afterafter
    derivative.push_back((-3.*data[size_t - 1][0] + 4.*data[size_t - 1][1] - data[size_t][2]) / (2. * space_step));

    //centered
    for (int i = 1; i < size_x - 1; i++)
        derivative.push_back((data[size_t - 1][i - 1] - data[size_t - 1][i + 1]) / (2. * space_step));

    //backward: 3*now - 4*before + beforebefore
    derivative.push_back((3*data[size_t - 1][size_x - 1] - 4*data[size_t - 1][size_x-2] + data[size_t-1][size_x-3]) / (2. * space_step));

    return derivative;
}

std::vector<double> WaveSolver::SecondDerSpaceCenteredDiff2(double **data, int size_t, int size_x, double space_step){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    std::vector<double> derivative;
    
    if(size_x < 4){
        std::cout << "[WARNING] The size of your data is not enough to compute this approximation, please try another one!" << std::endl;
        return derivative;
    }

    //forward
    derivative.push_back((2.*data[size_t - 1][0] - 5. * data[size_t - 1][1] + 4.* data[size_t - 1][2] - data[size_t-1][3]) / (space_step * space_step * space_step));

    //centered
    for (int i = 1; i < size_x - 1; i++)
        derivative.push_back((data[size_t - 1][i - 1] - 2. * data[size_t - 1][i] + data[size_t - 1][i + 1]) / (space_step * space_step));

    //backward
    derivative.push_back((2.*data[size_t - 1][size_x - 1] - 5. * data[size_t - 1][size_x - 2] + 4.*data[size_t - 1][size_x-3] - data[size_t-1][size_x-4]) / ( space_step * space_step * space_step));

    return derivative;
}

////////////////////////
////////////////////////
////////////////////////
//////WRITING DATA//////
////////////////////////
////////////////////////
////////////////////////

void WaveSolver::Write(std::string filename, double **data, int size_t, int size_x, double space_step, double time_step){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    
    std::ofstream myfile;
    myfile.open(filename.c_str());

    myfile << size_t << " " << size_x << " " << space_step << " " << time_step << "\n";

    for(int i = 0; i < size_t; i++){
        for(int j = 0; j < size_x; j++)
            myfile << data[i][j] << " ";
        myfile << "\n";
    }

    myfile.close();
}