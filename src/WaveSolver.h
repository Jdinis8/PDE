#ifndef __WaveSolver__
#define __WaveSolver__

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>


//#define DEBUG

class WaveSolver{
    public:
        WaveSolver();
        ~WaveSolver() = default;

        //In the following methods, we assume data[i][j] means ith step of time, jth step of space
        
        // Time Methods
        void TimeSimpleDiff1  (double**& idata, int size_t, int size_x, double space_step, double time_step, std::vector<double> RHS); 
        void TimeCenteredDiff2(double**& idata, int size_t, int size_x, double space_step, double time_step, std::vector<double> RHS);
        void TimeWaveRK4      (double **&u_data, double **&udot_data, int size_t, int size_x, double space_step, double time_step); //this solves the all wave equation immediately

        // Periodic Spatial Methods
        std::vector<double> PFirstDerSpaceCenteredDiff2 (double **idata, int size_t, int size_x, double space_step);
        std::vector<double> PFirstDerSpaceCenteredDiff4 (double **data, int size_t, int size_x, double space_step);
        std::vector<double> PSecondDerSpaceCenteredDiff2(double **data, int size_t, int size_x, double space_step);
        

        //Non-Periodic Spatial Methods
        std::vector<double> FirstDerSpaceCenteredDiff2 (double **data, int size_t, int size_x, double space_step);
        std::vector<double> SecondDerSpaceCenteredDiff2(double **data, int size_t, int size_x, double space_step);

        //Spectral Method
        std::vector<double> PseudoSpectral(double **data, int size_t, int size_x, double space_step);

        //Convergence
        std::vector<std::vector<double>> ConvergenceTest(double **udata, double **udotdata, int size_t, int size_x, double space_step, double time_step, int f, int order);
        std::vector<double>              L2NormTime(double **udata, double **udotdata, int size_t, int size_x, double space_step, double time_step, int f);
        double                           L2NormStep(double **udata, double **udotdata, int size_t, int size_x, double space_step, double time_step, int f, int order);


        // Writing stuff to a txt file
        void Write(std::string filename, double **data, int size_t, int size_x, double space_step, double time_step, double x0);
        void Write(std::string filename, std::vector<std::vector<double>> data, int size_t, int size_x, double space_step, double time_step, double x0); //just to write the convergence data
        void Write(std::string filename, double **data, int size_h);
        void Write(std::string filename, std::vector<double> data, double time_step);
};

#endif