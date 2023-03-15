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

    for (int i = 0; i < size_t; i++) delete[] dummy[i];
    delete[] dummy;
}

void WaveSolver::TimeWaveRK4(double**& u_data, double**& udot_data, int size_t, int size_x, double space_step, double time_step){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    //this one must be a more specific method, so we will do this in a not so general way, but easy to adapt, in case it's necessary
    //we'll have to have two arrays, one for the time steps with h and another
    //one for time steps of h/2

    //furthermore, we will assume we have a system as:
    //dot(pi) = u''
    //dot(u) = pi
    //which is equal to: u'' = dot(dot(u))

    //we will also assume that the data is being given in the following way
    //first n places are for u(x,t) and the other n places are for pi(x,t)
    //i.e., normalstep[i][0-N] is u, normalstep[i][N-2] is for pi

    std::vector<std::vector<double>> K1, K2, K3, K4;
    std::vector<double> dummy1;
    
    double **new_udata    = new double*[size_t + 1];
    double **new_udotdata = new double*[size_t + 1];

    for (int j = 0; j < size_t + 1; j++)
    {
        new_udata[j] = new double[size_x];
        new_udotdata[j] = new double[size_x];
    }

    for (int k = 0; k < size_x; k++){
        new_udata[0][k] = u_data[0][k];
        new_udotdata[0][k] = udot_data[0][k];
    }

    for (int i = 0; i < size_t; i++){
        K1.clear();
        K2.clear();
        K3.clear();
        K4.clear();
        dummy1.clear();


        //K1[0] = u'', K1[1] = pi
        for (int j = 0; j < size_x; j++) dummy1.push_back(new_udotdata[i][j]);
        
        K1.push_back(PSecondDerSpaceCenteredDiff2(new_udata, i+1, size_x, space_step));
        K1.push_back(dummy1);

        
        //we briefly change udata for K2
        for(int j = 0; j < size_x; j++){
            new_udata[i][j] += time_step*K1[0][j]*0.5;
                  dummy1[j] += time_step*K1[1][j]*0.5;
        }

        K2.push_back(PSecondDerSpaceCenteredDiff2(new_udata, i+1, size_x, space_step));
        K2.push_back(dummy1);

        //we subtract K1 and add K2 to get u_data+K2
        for(int j = 0; j < size_x; j++){
            new_udata[i][j] += time_step*0.5*(K2[0][j]-K1[0][j]);
                  dummy1[j] += time_step*0.5*(K2[1][j]-K1[1][j]);
        }

        K3.push_back(PSecondDerSpaceCenteredDiff2(new_udata, i+1, size_x, space_step));
        K3.push_back(dummy1);

        //we do the same, but this time we take care to subtract K2/2 and add K3
        for(int j = 0; j < size_x; j++){
            new_udata[i][j] += time_step*(K3[0][j]-K2[0][j]*0.5);
                  dummy1[j] += time_step*(K3[1][j]-K2[1][j]*0.5);
        }

        K4.push_back(PSecondDerSpaceCenteredDiff2(new_udata, i+1, size_x, space_step));
        K4.push_back(dummy1);

        //correcting new_udata[i][j] to being as it was before
        for (int j = 0; j < size_x; j++) new_udata[i][j] += time_step * (-K3[0][j]);

        //now that we have calculated all k's, it is time to finally calculate the next step in time
        //recall: new_udotdata = pi,  udata = u; thus, pi_dot = u'', u_dot = pi;
        for(int j = 0; j < size_x; j++){
            new_udotdata[i+1][j] = new_udotdata[i][j] + time_step/6.*(K1[0][j] + 2.*K2[0][j] + 2.*K3[0][j] + K4[0][j]);
               new_udata[i+1][j] = new_udata[i][j] + time_step/6.*(K1[1][j] + 2.*K2[1][j] + 2.*K3[1][j] + K4[1][j]);
        }
    }

    double **dumb = u_data;
    double **dumb2 = udot_data;

    u_data = new_udata;
    udot_data = new_udotdata;

    delete[] dumb;
    //delete[] dumb2;
    
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

std::vector<double> WaveSolver::PseudoSpectral(double **data, int size_t, int size_x, double space_step){

}

////////////////////////
////////////////////////
////////////////////////
//////WRITING DATA//////
////////////////////////
////////////////////////
////////////////////////

void WaveSolver::Write(std::string filename, double **data, int size_t, int size_x, double space_step, double time_step, double x0){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    
    std::ofstream myfile;
    myfile.open(filename.c_str());

    myfile << size_t << " " << size_x << " " << space_step << " " << time_step << " " << x0 << "\n";

    for(int i = 0; i < size_t; i++){
        for(int j = 0; j < size_x; j++)
            myfile << data[i][j] << " ";
        myfile << "\n";
    }

    myfile.close();
}