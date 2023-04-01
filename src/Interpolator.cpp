#include "Interpolator.h"

//#define DEBUG

Interpolator::Interpolator(double *ix, double *iy, int iN): x(ix), y(iy), N(iN){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
}

////////////////////
////////////////////
//////LAGRANGE//////
////////////////////
////////////////////

double Interpolator::Lagrange(double fx){
	#ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

	double fy=0;
	for(int i=0; i<N; i++){
		double l = 1;
		for(int j=0; j<N; j++){
			if(j != i)
				l *= (fx - x[j]) / (x[i] - x[j]);
		}

		fy += l * y[i];
	}

	return fy;
}

double Interpolator::LagrangeFirstDerivative(double fx){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    if(N==1) return 0;
	if(fx == x[0] || fx == x[N-1]){
        std::cout << "[WARNING] You can't calculate FirstDerivative on the first point nor the last!" << std::endl;
        return 0.;
    }

    return (Lagrange(fx + h) - Lagrange(fx - h)) / (2 * h);
}

double Interpolator::LagrangeSecondDerivative(double fx){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    if(N==1) return 0;
	if(fx == x[0] || fx == x[N-1]){
        std::cout << "[WARNING] You can't calculate FirstDerivative on the first point nor the last!" << std::endl;
        return 0.;
    }

    return (Lagrange(fx + h) - 2. * Lagrange(fx) + Lagrange(fx - h)) / (h * h);
}

////////////////////
////////////////////
///////NEWTON///////
////////////////////
////////////////////

double Interpolator::Newton(double fx){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
    //setting up coefficientes first
    double* coef = new double[N];
    for(int i=0; i<N; i++) coef[i] = y[i];
	for(int i=1; i<N; i++){
		for(int j=i; j<N; j++)
			coef[j] = (coef[j] - coef[i-1]) / (x[j] - x[i-1]);
	}
	
    //interpolating
	double fy=coef[N-1];
	for(int i=1; i<N; i++) fy = coef[N-1-i] + (fx - x[N-1-i]) * fy;

    delete[] coef;
    return fy;
}

double Interpolator::NewtonFirstDerivative(double fx){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

    if(N==1) return 0;
	if(fx == x[0] || fx == x[N-1]){
        std::cout << "[WARNING] You can't calculate SecondDerivative on the first point nor the last!" << std::endl;
        return 0.;
    }

    return (Newton(fx + h) - Newton(fx - h)) / (2 * h);
}

double Interpolator::NewtonSecondDerivative(double fx){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

    if(N==1) return 0;
	if(fx == x[0] || fx == x[N-1]){
        std::cout << "[WARNING] You can't calculate SecondDerivative on the first point nor the last!" << std::endl;
        return 0.;
    }

    return (Newton(fx + h) - 2. * Newton(fx) + Newton(fx - h)) / (h * h);
}

////////////////////
////////////////////
//////NEVILLE///////
////////////////////
////////////////////

double Interpolator::Neville(double fx){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

	double* p = new double[N];
	for(int i=0; i<N; i++)
		p[i] = y[i];

	for(int k=1; k<N; k++){
		for(int i=0; i<N-k; i++)
			p[i] = ((fx - x[i+k]) * p[i] - (fx - x[i]) * p[i+1]) / (x[i] - x[i+k]);
	}
	
	double fy = p[0];
	delete p;

    return fy;
}

double Interpolator::NevilleFirstDerivative(double fx){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

    if(N==1) return 0;
	if(fx == x[0] || fx == x[N-1]){
        std::cout << "[WARNING] You can't calculate FirstDerivative on the first point nor the last!" << std::endl;
        return 0.;
    }

    return (Neville(fx + h) - Neville(fx - h)) / (2 * h);
}

double Interpolator::NevilleSecondDerivative(double fx){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

    if(N==1) return 0;
	if(fx == x[0] || fx == x[N-1]){
        std::cout << "[WARNING] You can't calculate SecondDerivative on the first point nor the last!" << std::endl;
        return 0.;
    }

    return (Neville(fx + h) - 2. * Neville(fx) + Neville(fx - h)) / (h * h);
}

void Interpolator::SetX(double* ix){
	#ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
	for(int i = 0; i < N; i++) this->x[i] = ix[i];
}

void Interpolator::SetY(double* iy){
	#ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
	for(int i = 0; i < N; i++) this->y[i] = iy[i];
}