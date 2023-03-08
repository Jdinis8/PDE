#ifndef __WaveSolver__
#define __WaveSolver__

#include <iostream>
#include <vector>

class WaveSolver{
    public:
        WaveSolver();
        ~WaveSolver() = default;

        //In the following methods, we assume data[i][j] means ith step of time, jth step of space
        
        // Time Methods
        void TimeSimpleDiff1  (double** data, int size_t, int size_x, double space_step, double time_step, std::vector<double> RHS); 
        void TimeCenteredDiff2(double **data, int size_t, int size_x, double space_step, double time_step, std::vector<double> RHS);

        //Periodic Spatial Methods
        std::vector<double> PFirstDerSpaceCenteredDiff2 (double **data, int size_t, int size_x, double space_step);
        std::vector<double> PFirstDerSpaceCenteredDiff4 (double **data, int size_t, int size_x, double space_step);
        std::vector<double> PSecondDerSpaceCenteredDiff2(double **data, int size_t, int size_x, double space_step);
        

        //Non-Periodic Spatial Methods
        std::vector<double> FirstDerSpaceCenteredDiff2 (double **data, int size_t, int size_x, double space_step);
        std::vector<double> SecondDerSpaceCenteredDiff2(double **data, int size_t, int size_x, double space_step);

        std::vector<double> RK4(double **data, int size_t, int size_x, double space_step, double time_step);
};

#endif