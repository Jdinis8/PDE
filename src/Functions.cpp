#include "Functions.h"

double f(double x, double t){
    return 0.;
} //this is the source function in the wave equation

double initialf(double x){
    return exp(-x*x/0.1);
}

double initialdf(double x){
    return 0;
}