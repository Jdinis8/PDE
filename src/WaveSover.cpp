#include "WaveSolver.h"

WaveSolver::WaveSolver(double* input_data, int isize_data, int it0, int itf, double speed): size_data(isize_data), t0(it0), tf(itf), c(speed){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    for (int i = 0; i < size_data; i++) initial_data.push_back(input_data[i]); 
}

std::vector<std::vector<double>> WaveSolver::CenteredDiff(double space_step, double time_step){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    double time_counter(t0);
    int    j(0); //j is time, i is space

    std::vector<std::vector<double>> Sol;
    Sol.push_back(initial_data);
    
    while(time_counter <= tf){
        std::vector<double> next;

        //for when we only have one time in the vector
        if(j == 0){

            j++;
        }else{

            //for when we are at the edges of the vector
            

            //for when we are in the middle
            for(int i = 1; i < initial_data.size()-1; i++){
                next.push_back(2 * Sol[i][j] * (1 - c * c * time_step * time_step / (space_step * space_step)) - Sol[i][j - 1] + c * c * time_step * time_step / (space_step * space_step)*(Sol[i+1][j] + Sol[i-1][j]));
            }
        }


        time_counter += time_step;
        j++;
    }

    return Sol;
}