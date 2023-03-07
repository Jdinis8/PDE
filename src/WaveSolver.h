#ifndef __WaveSolver__
#define __WaveSolver__

#include <iostream>
#include <vector>

class WaveSolver{
    public:
        WaveSolver(double& input_data, int size_data, int it0, int itf);
        ~WaveSolver() = default;

        std::vector<std::vector<double>> CenteredDiff(double space_step, double time_step);
        std::vector<std::vector<double>> BackwardDiff(double space_step, double time_step);
        std::vector<std::vector<double>> ForwardDiff (double space_step, double time_step);

    private:
        double initial_data;
        int size_data, t0, tf;
};

#endif