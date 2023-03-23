#include "Matrix.h"

Matrix::Matrix(double** imtx, int im, int in): mtx(imtx), m(im), n(in){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
}

Matrix::Matrix(std::vector<std::vector<double>> imtx): m(int(imtx.size())), n(int(imtx[0].size())){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    mtx = new double*[m];
    for(int i = 0; i < m; i++){
        mtx[i] = new double[n];
        for(int j = 0; j < n; j++) mtx[i][j] = imtx[i][j];
    }
}

Matrix::~Matrix(){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    for(int i = 0; i < m; i++) delete[] mtx[i];
    delete[] mtx;
}

double* Matrix::Mult(double* vec){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    double* res = new double[m];
    for(int i = 0; i < m; i++){
        res[i] = 0.;
        for(int j = 0; j < n; j++) res[i] += mtx[i][j]*vec[j];
    }
    return res;
}

std::vector<double> Matrix::Mult(std::vector<double> vec){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    std::vector<double> res;
    double sum(0.);
    for(int i = 0; i < m; i++){
        sum = 0.;
        for(int j = 0; j < n; j++) sum += mtx[i][j]*vec[j];
        res.push_back(sum);
    }
    return res;
}

std::vector<std::vector<double>> Matrix::Mult(std::vector<std::vector<double>> mtx2){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    std::vector<std::vector<double>> res;
    std::vector<double> dumb;
    double sum(0.);

    for(int i = 0; i < m; i++){
        dumb.clear();
        for(int j = 0; j < int(mtx2.size()); j++){
            sum = 0.;
            for(int k = 0; k < n; k++) sum += mtx[i][k]*mtx2[k][j];
            dumb.push_back(sum);
        }
        res.push_back(dumb);
    }
    return res;
}

Matrix* Matrix::ChebyshevDN(double *chebyshevpts, int size){
    double** mtp = new double*[size];
    for (int i = 0; i < size; i++) mtp[i] = new double[size];

    
    Matrix* mtx = new Matrix(mtp, size, size);

    return mtx;
}

    /////////////////
    /////////////////
    ////OPERATORS////
    /////////////////
    /////////////////

    std::vector<double> &Matrix::operator[](int i)
{
#ifdef DEBUG
    printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
    std::vector<double>* vec = new std::vector<double>;
    if(i < 0 || i >= m){
        std::cout << "[WARNING] You are trying to access a row of the matrix which does not existe!!" << std::endl;
        return *vec;
    }
    for(int j = 0; j < m; j++) vec->push_back(mtx[i][j]);
    return *vec;
}

std::vector<double> Matrix::operator[](int i) const{
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    std::vector<double> vec;
    if(i < 0 || i >= m){
        std::cout << "[WARNING] You are trying to access a row of the matrix which does not exist!!" << std::endl;
        return vec;
    }
    for(int j = 0; j < m; j++) vec.push_back(mtx[i][j]);
    return vec;
}

int Matrix::GetRowN() const{return m;}
int Matrix::GetColN() const{return n;}

void Matrix::Print() const{
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    std::cout << "Matrix: " << std::endl;
	for (int i = 0; i < m; i++){
        std::cout << "              ";
        for(int j = 0; j < n; j++)
		    std::cout << mtx[i][j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

std::ostream& operator<<(std::ostream& s, const Matrix& mtx){
	s << "Matrix Full: " << std::endl;
  	for (int i=0; i< mtx.GetRowN() ; i++){
  		s << "              ";
        for(int j = 0; j < mtx.GetColN(); j++)
            s << mtx[i][j] << " ";
        s << std::endl;
    }
 	 return s;
}