#ifndef __Functions__
#define __Functions__

#include <math.h>

//Source Function for the Wave Equation
double f(double x, double t);

//Initial data for the wave equation: f(x,0) and d_t(f)(x,0)
double initialf (double x);
double initialdf(double x);

#endif