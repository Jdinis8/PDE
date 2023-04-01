#ifndef __Interpolator__
#define __Interpolator__

#include <iostream>

//#define DEBUG

class Interpolator{
    public:
        Interpolator(double *ix, double *iy, int iN);
        ~Interpolator() = default;
        
        //Lagrange Interpolator Methods
        double Lagrange(double fx);
        double LagrangeFirstDerivative(double fx);
        double LagrangeSecondDerivative(double fx);

        //Newton Interpolator Methods
        double Newton(double fx);
        double NewtonFirstDerivative(double fx);
        double NewtonSecondDerivative(double fx);

        //Neville Interpolator Methods
        double Neville(double fx);
        double NevilleFirstDerivative(double fx);
        double NevilleSecondDerivative(double fx);

        //Change Data for Interpolation
        void SetX(double* ix);
        void SetY(double* iy);

    private:
        double* x;
        double* y;
        int N;
        double h;
};

#endif