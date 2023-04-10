#ifndef __WaveSolver__
#define __WaveSolver__

#include <sstream>
#include <fstream>
#include "Matrix.h"
#include "Functions.h"

class WaveSolver{
    public:
        WaveSolver(double, int);
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
        std::vector<double> PseudoSpectral   (double **data, int size_t, int size_x);
        std::vector<double> BVPseudoSpectral (double** data, int size_t, int size_x, double initial, double final);
        std::vector<double> FFTPseudoSpectral(double** data, int size_t, int size_x);
        std::vector<double> FFTPseudoSpectral(std::vector<double> data, int size_x);

        //Convergence
        std::vector<std::vector<double>> PointConvergenceTest(double **udata, double **udotdata, int size_t, int size_x, double space_step, double time_step, int f, int order);
        std::vector<double>              L2NormTime(double **udata, double **udotdata, int size_t, int size_x, double space_step, double time_step, int f);
        double                           L2NormStep(double **udata, double **udotdata, int size_t, int size_x, double space_step, double time_step, int f, int order);

        //ghost propagator
        void Ghost(double** udata, double** udotdata, int& size_t, int& size_x); //adds amt_ghost points to both sides of the last element of the array udata and udotdata

        //method selector
        void        SetMethod(std::string imethod);
        std::string GetMethod();

        // Writing stuff to a txt file
        void Write(std::string filename, double **data, int size_t, int size_x, double space_step, double time_step, double x0);
        void Write(std::string filename, double **data, int size_t, int size_x, double time_step, double* x);
        void Write(std::string filename, std::vector<std::vector<double>> data, int size_t, int size_x, double space_step, double time_step, double x0); //just to write the convergence data
        void Write(std::string filename, double **data, int size_h);
        void Write(std::string filename, std::vector<double> data, double time_step);
    
    private:
        double CFL;
        int amt_ghost;
        std::string method;
};

#endif