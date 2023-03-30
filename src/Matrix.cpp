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

Matrix::Matrix(): m(0), n(0){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
}

std::vector<double> Matrix::Mult(double* vec){
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

Matrix Matrix::ChebyshevDN(int N){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif

    double* chebyshevpts = new double[N];
    for(int i = 0; i < N; i++) chebyshevpts[i] = cos(i*M_PI/N);

    double** mtp = new double*[N+1];
    for (int i = 0; i < N+1; i++) mtp[i] = new double[N+1];
    double ci(1.), cj(1.);

    //we first populate the four edges of the matrix
    mtp[0][0] = (2*N*N+1.)/6.;
    mtp[N][N] = -mtp[0][0];
    mtp[N][0] = -0.5*pow(-1, N);
    mtp[0][N] = -mtp[N][0];

    //then the diagonal
    for(int i = 1; i < N; i++) mtp[i][i]= -chebyshevpts[i]/(2.*(1.-chebyshevpts[i]*chebyshevpts[i]));

    //the off terms of the diagonal, except the N+1'th row and the N+1'th column
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(i != j){
                if(i == 0 || i == N) ci = 2.;
                else ci = 1.;
                if(j == 0 || j == N) cj = 2.;
                else cj = 1.;
                mtp[i][j] = ci/cj*pow(-1, i + j) / (chebyshevpts[i] - chebyshevpts[j]);
            }
        }
    }

    for(int i = 1; i < N; i++){
        mtp[i][N] = 0.5*pow(-1, N+i)/(1.+chebyshevpts[i]);
        mtp[N][i] =-2.*pow(-1, N+i)/(1.+chebyshevpts[i]);
    }

    delete[] chebyshevpts;
    return Matrix(mtp, N+1, N+1);
}

/////////////////
/////////////////
////OPERATORS////
/////////////////
/////////////////

std::vector<double> &Matrix::operator[](int i){
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
    //#ifdef DEBUG
    //    printf("[%s]\n", __PRETTY_FUNCTION__);
    //#endif
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

Matrix Matrix::operator*(const Matrix& mtx2){
    #ifdef DEBUG
        printf("[%s]\n", __PRETTY_FUNCTION__);
    #endif
    std::vector<std::vector<double>> res;
    std::vector<double> dumb;
    double sum(0.);

    for(int i = 0; i < m; i++){
        dumb.clear();
        for(int j = 0; j < int(mtx2.GetColN()); j++){
            sum = 0.;
            for(int k = 0; k < n; k++) sum += mtx[i][k]*mtx2[k][j];
            dumb.push_back(sum);
        }
        res.push_back(dumb);
    }
    return Matrix(res);
}