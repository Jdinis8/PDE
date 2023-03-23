#ifndef __Matrix__
#define __Matrix__

#include <vector>
#include <iostream>
#include <math.h>

//#define DEBUG

class Matrix{
    public:
        Matrix(double** imtx, int im, int in);
        Matrix(std::vector<std::vector<double>> imtx);
        ~Matrix();

        double*                          Mult(double* vec);
        std::vector<double>              Mult(std::vector<double> vec);
        double**                         Mult(double** mtx2, int m2, int n2);
        std::vector<std::vector<double>> Mult(std::vector<std::vector<double>> mtx2);

        Matrix ChebyshevDN(double* chebyshevpts, int N);

        //to get info from the matrix
        int GetRowN() const;
        int GetColN() const;

        //operators and printers
        std::vector<double> &operator[](int i);
        std::vector<double>  operator[](int i) const;

        void Print() const;
        friend std::ostream &operator<<(std::ostream &, const Matrix &);

    private: 
        double **mtx;
        int m, n; //m lines, n columns
};

#endif