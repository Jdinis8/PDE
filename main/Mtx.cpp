#include <iostream>
#include "Matrix.h"

int main(){

    std::vector<std::vector<double>> imtx;
    std::vector<double> line1;
    std::vector<double> line2;

    line1.push_back(1);
    line1.push_back(2);

    line2.push_back(2);
    line2.push_back(1);

    imtx.push_back(line1);
    imtx.push_back(line2);

    std::vector<std::vector<double>> imtx2;
    std::vector<double> line12;
    std::vector<double> line22;

    line12.push_back(1);
    line12.push_back(0);

    line22.push_back(0);
    line22.push_back(1);

    imtx2.push_back(line12);
    imtx2.push_back(line22);

    std::vector<double> vec;
    vec.push_back(0);
    vec.push_back(1);

    Matrix m1(imtx);

    std::vector<std::vector<double>> mmult = m1.Mult(imtx2);

    
    std::cout << mmult[0][0] << " " << mmult[0][1] << std::endl;
    std::cout << mmult[1][0] << " " << mmult[1][1] << std::endl;
    
    std::vector<double> vmult = m1.Mult(vec);

    std::cout << vmult[0] << " " << vmult[1] << std::endl;

    int ncheb(5);
    double* cheb = new double[ncheb];
    for(int i = 0; i < ncheb; i++) cheb[i] = cos(i*M_PI/ncheb);
    
    Matrix Mcheb = m1.ChebyshevDN(cheb, ncheb);

    std::cout << Mcheb << std::endl;


    return 0;
}