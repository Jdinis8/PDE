#include "WaveSolver.h"

WaveSolver::WaveSolver(double &input_data, int isize_data, int it0, int itf): size_data(isize_data), t0(it0), tf(itf), initial_data(input_data){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
}

std::vector<std::vector<double>> WaveSolver::CenteredDiff(){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

}