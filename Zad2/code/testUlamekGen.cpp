#include <stdio.h>
#include <iostream>
#include <list>
#include <vector>
#include <iomanip>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <thread>
#include <Eigen/Dense>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include "gauss.h"
#define MAXU 21
#define MAX 801
#define PRECD 16
#define PRECF 8
using namespace Eigen;
using namespace std;
//g++ -O3 -I eigen-eigen-5a0156e40feb -std=c++11 testUlamekGen.cpp -o testu.exe
int NWD(int a, int b)
{
    if(a<0) cout<< "no chyba nie\n";
    if(b<0) cout<< "ty no nie wiem\n";
    while(a!=b)
       if(a>b)
           a-=b;
       else
           b-=a;
    return a;
}

unsigned long long gcd(unsigned long long a, unsigned long long b)
{
    for (;;)
    {
        if (a == 0) return b;
        b %= a;
        if (b == 0) return a;
        a %= b;
    }
}

unsigned long long lcm(unsigned long long a, unsigned long long b)
{
    unsigned long long temp = gcd(a, b);

    return temp ? (a / temp * b) : 0;
}

class Fraction {
    public:
        Fraction(unsigned long long num, unsigned long long denum) : _num(num), _denum(denum) {}
        Fraction(){}

        void set_num(unsigned long long x){
            _num=x;
        }
        unsigned long long get_num(){
            return _num;
        }
        void set_denum(unsigned long long x){
            _denum=x;
        }
        unsigned long long get_denum(){
            return _denum;
        }

        void operator*= (const Fraction& other) {
            this->_num *= other._num;
            this->_denum *= other._denum;
        }

        void operator/= (const Fraction& other) {
            this->_num *= other._denum;
            this->_denum *= other._num;
        }

        void operator+= (const Fraction& other) {
            unsigned long long arr[] = { this->_denum, other._denum};
            unsigned long long  res = std::accumulate(arr, arr + 2, 1, lcm);
            this->_num = (this->_num * res/this->_denum) + other._num * res/other._denum;
            this->_denum = res;
        }

        void operator-= (const Fraction& other) {
            unsigned long long arr[] = { this->_denum, other._denum};
            unsigned long long  res = std::accumulate(arr, arr + 2, 1, lcm);
            this->_num = (this->_num * res/this->_denum) - other._num * res/other._denum;
            this->_denum = res;
        }
    private:
        unsigned long long _num;
        unsigned long long _denum;
};

Fraction operator* (const Fraction& lhs, const Fraction& rhs) {
    Fraction temp (lhs);
    temp *= rhs;
    int x = NWD(temp.get_num(), temp.get_denum());
    temp.set_num(temp.get_num()/x);
    temp.set_denum(temp.get_denum()/x);
    return temp;
}

Fraction operator/ (const Fraction& lhs, const Fraction& rhs) {
    Fraction temp (lhs);
    temp /= rhs;
    return temp;
}

Fraction operator+ (const Fraction& lhs, const Fraction& rhs) {
    Fraction temp (lhs);
    temp += rhs;
    return temp;
}

Fraction operator- (const Fraction& lhs, const Fraction& rhs) {
    Fraction temp (lhs);
    temp -= rhs;
    return temp;
}


template <typename T> class MyMatrix {
public:
    vector<vector<T> > matrix;
private:
    unsigned rows;
    unsigned cols;

    T abs(T x) {
        return (x >= 0) ? x : -x;
    }

    public:
        // konstruktory
        MyMatrix(unsigned wiersz, unsigned kolumna, const T& x) {
            matrix.resize(wiersz);
            for (unsigned i=0; i<matrix.size(); i++) {
                matrix[i].resize(kolumna, x);
            }
            rows = wiersz;
            cols = kolumna;
        }
        MyMatrix(const MyMatrix<T>& pom) {
            matrix = pom.matrix;
            rows = pom.getRowCount();
            cols = pom.getColCount();
        }


        MyMatrix<T>& operator=(const MyMatrix<T>& pom) {
            if (&pom == this) return *this;
            unsigned new_rows = pom.getRowCount();
            unsigned new_cols = pom.getColCount();

            matrix.resize(new_rows);
            for (unsigned i=0; i<matrix.size(); i++) {
                matrix[i].resize(new_cols);
            }

            for (unsigned i=0; i<new_rows; i++) {
                for (unsigned j=0; j<new_cols; j++) {
                    matrix[i][j] = pom(i, j);
                }
            }
            rows = new_rows;
            cols = new_cols;
            return *this;
        }

        MyMatrix<T> operator+(const MyMatrix<T>& pom) {
            MyMatrix result(rows, cols);
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    result(i,j) = this->matrix[i][j] + pom(i,j);
                }
            }
            return result;
        }
        MyMatrix<T>& operator+=(const MyMatrix<T>& pom) {
            unsigned rows = pom.getRowCount();
            unsigned cols = pom.getColCount();
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    this->matrix[i][j] += pom(i,j);
                }
            }
            return *this;
        }
        MyMatrix<T> operator-(const MyMatrix<T>& pom) {
            MyMatrix result(rows, cols);
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    result(i,j) = this->matrix[i][j] - pom(i,j);
                }
            }
            return result;
        }
        MyMatrix<T>& operator-=(const MyMatrix<T>& pom) {
            unsigned rows = pom.getRowCount();
            unsigned cols = pom.getColCount();
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    this->matrix[i][j] -= pom(i,j);
                }
            }
            return *this;
        }
        MyMatrix<T> operator*(const MyMatrix<T>& pom) {
            unsigned rows = pom.getRowCount();
            unsigned cols = pom.getColCount();
            MyMatrix result(rows, cols);
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    for (unsigned k=0; k<rows; k++) {
                        result(i,j) += this->matrix[i][k] * pom(k,j);
                    }
                }
            }
            return result;
        }
        MyMatrix<T>& operator*=(const MyMatrix<T>& pom) {
            MyMatrix result = (*this) * pom;
            (*this) = result;
            return *this;
        }
        MyMatrix<T> transpose() {
            MyMatrix result(rows, cols);
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    result(i,j) = this->matrix[j][i];
                }
            }
            return result;
        }

        MyMatrix<T> operator+(const T& pom) {
            MyMatrix result(rows, cols);
            for (unsigned i = 0; i < rows; i++) {
                for (unsigned j = 0; j < cols; j++) {
                    result(i,j) = this->matrix[i][j] + pom;
                }
            }
            return result;
        }
        MyMatrix<T> operator-(const T& pom) {
            MyMatrix result(rows, cols);
            for (unsigned i = 0; i < rows; i++) {
                for (unsigned j = 0; j < cols; j++) {
                    result(i,j) = this->matrix[i][j] - pom;
                }
            }
            return result;
        }
        MyMatrix<T> operator*(const T& pom) {
            MyMatrix result(rows, cols);
            for (unsigned i = 0; i < rows; i++) {
                for (unsigned j = 0; j < cols; j++) {
                    result(i,j) = this->matrix[i][j] * pom;
                }
            }
            return result;
        }
        MyMatrix<T> operator/(const T& pom) {
            MyMatrix result(rows, cols);
            for (unsigned i = 0; i < rows; i++) {
                for (unsigned j = 0; j < cols; j++) {
                    result(i,j) = this->matrix[i][j] / pom;
                }
            }
            return result;
        }

        vector<T> operator*(const vector<T>& pom) {
            vector<T> result(pom.size());
            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    result[i] += this->matrix[i][j] * pom[j];
                }
            }
            return result;
        }
        vector<T> diagonalVector() {
            vector<T> result(rows);
            for (unsigned i=0; i<rows; i++) {
                result[i] = this->matrix[i][i];
            }
            return result;
        }

        T& operator()(const unsigned& row, const unsigned& col) {
            return this->matrix[row][col];
        }
        const T& operator()(const unsigned& row, const unsigned& col) const {
            return this->matrix[row][col];
        }

        void setAt(unsigned row, unsigned col, const T& x) {
            this->matrix[row][col] = x;
        }

        T& getAt(unsigned row, unsigned col) {
            return this->matrix[row][col];
        }

        unsigned getRowCount() const {
            return this->rows;
        }
        unsigned getColCount() const {
            return this->cols;
        }

//Gauss pełnego wyboru
        vector<T> Gauss1() {
            int n = getRowCount();
            for (int i = 0; i < n; i++) {
                T maxEl = abs(matrix[i][i]);
                int maxRow = i;
                int maxCol = i;
                for (int k = i+1; k < n; k++) {
                    if (abs(matrix[k][i]) > maxEl) {
                        maxEl = abs(matrix[k][i]);
                        maxRow = k;
                        maxCol = i;
                    }
                }
                for (int k = i; k < n+1; k++) {
                    T pom = matrix[maxRow][k];
                    matrix[maxRow][k] = matrix[i][k];
                    matrix[i][k] = pom;
                }
                for (int k = i; k < n; k++) {
                    T pom = matrix[k][maxCol];
                    matrix[k][maxCol] = matrix[k][i];
                    matrix[k][i] = pom;
                }

                for (int k = i+1; k < n; k++) {
                    T c = -matrix[k][i] / matrix[i][i];
                    for (int j = i; j < n+1; j++) {
                        if (i == j) {
                            matrix[k][j] = 0;
                        } else {
                            matrix[k][j] += c * matrix[i][j];
                        }
                    }
                }
            }
            vector<T> x(n);
            for (int i=n-1; i>=0; i--) {
                x[i] = matrix[i][n] / matrix[i][i];
                for (int k=i-1;k>=0; k--) {
                    matrix[k][n] -= matrix[k][i] * x[i];
                }
            }
            return x;
        }



   //Gauss połowiczny wybr
        vector<T> Gauss2() {
            int n = getRowCount();
            for (int i = 0; i < n; i++) {
                T maxEl = abs(matrix[i][i]);
                int maxRow = i;
                for (int k = i+1; k < n; k++) {
                    if (abs(matrix[k][i]) > maxEl) {
                        maxEl = abs(matrix[k][i]);
                        maxRow = k;
                    }
                }
                for (int k = i; k < n+1; k++) {
                    T pom = matrix[maxRow][k];
                    matrix[maxRow][k] = matrix[i][k];
                    matrix[i][k] = pom;
                }
                for (int k = i+1; k < n; k++) {
                    T c = -matrix[k][i] / matrix[i][i];
                    for (int j = i; j < n+1; j++) {
                        if (i == j) {
                            matrix[k][j] = 0;
                        } else {
                            matrix[k][j] += c * matrix[i][j];
                        }
                    }
                }
            }

            vector<T> x(n);
            for (int i=n-1; i>=0; i--) {
                x[i] = matrix[i][n] / matrix[i][i];
                for (int k=i-1;k>=0; k--) {
                    matrix[k][n] -= matrix[k][i] * x[i];
                }
            }
            return x;
        }

        //bez wyboru
          vector<T> Gauss3() {
            int n = getRowCount();
            for (int i = 0; i < n; i++) {


                // wyprowadź zera przed obecnym wierszem
                for (int k = i+1; k < n; k++) {
                    T c = -matrix[k][i] / matrix[i][i];
                    for (int j = i; j < n+1; j++) {
                        if (i == j) {
                            matrix[k][j] = 0;
                        } else {
                            matrix[k][j] += c * matrix[i][j];
                        }
                    }
                }
            }
            // rozwiąż Ax = B za pomocą powstałej macierzy trójkątnej
            vector<T> x(n);
            for (int i=n-1; i>=0; i--) {
                x[i] = matrix[i][n] / matrix[i][i];
                for (int k=i-1;k>=0; k--) {
                    matrix[k][n] -= matrix[k][i] * x[i];
                }
            }
            return x;
        }
};

void ulamek(){
    int i, t, x, y, h=0, m;  // liczniki
    std::ifstream fp("macierze.txt");
    ofstream abcXCzas;
    abcXCzas.open ("FactorCzasyABCX.txt");
    abcXCzas.precision(PRECD);

    ofstream aXCzas;
    aXCzas.open ("FactorCzasyAX.txt");
    aXCzas.precision(PRECD);

    ofstream abcCzas;
    abcCzas.open ("FactorCzasyABC.txt");
    abcCzas.precision(PRECD);

    ofstream abcXWyn;
    abcXWyn.open ("FactorWynikiABCX.txt");
    abcXWyn.precision(PRECD);

    ofstream aXWyn;
    aXWyn.open ("FactorWynAX.txt");
    aXWyn.precision(PRECD);

    ofstream abcWyn;
    abcWyn.open ("FactorWynABC.txt");
    abcWyn.precision(PRECD);


    for(i=10;i<MAXU;i = i+ 10){
        Fraction f2dArrayA[i][i];
        Fraction f2dArrayB[i][i];
        Fraction f2dArrayC[i][i];
        Fraction fVector[i][1];
        Fraction fWyn[i][i];
        Fraction fWyn2[i][i];
        Fraction fWyn3[i][i];
        Fraction fWynVector[i][1];
        Fraction fWynVector2[i][1];
        Fraction fWynAbc[i][i];

        for ( x=0; x<i;x++){
            fWynVector[x][0].set_num(0);
            fWynVector[x][0].set_denum(1);
            fWynVector2[x][0].set_num(0);
            fWynVector2[x][0].set_denum(1);
            for(y=0; y<i;y++){
                fWyn3[x][y].set_num(0);
                fWyn3[x][y].set_denum(1);
                fWyn2[x][y].set_num(0);
                fWyn2[x][y].set_denum(1);
                fWyn[x][y].set_num(0);
                fWyn[x][y].set_denum(1);
                fWynAbc[x][y].set_num(0);
                fWynAbc[x][y].set_denum(1);
                unsigned  long long r1;
                fp >> r1;
                unsigned  long long r2;
                fp >> r2;
                f2dArrayA[x][y].set_num(r1);
                f2dArrayA[x][y].set_denum(r2);
                fp >> r1;
                fp >> r2;
                f2dArrayB[x][y].set_num(r1);
                f2dArrayB[x][y].set_denum(r2);
                fp >> r1;
                fp >> r2;
                f2dArrayC[x][y].set_num(r1);
                f2dArrayC[x][y].set_denum(r2);
            }
        }

        for (t=0; t<i; t++) { // wektor
            unsigned long long r1;
            fp >> r1;
            unsigned long long r2;
            fp >> r2;
            fVector[t][0].set_num(r1);
            fVector[t][0].set_denum(r2);
        }

///////////////////////////(A + B + C) * X//////////////////////
        abcXCzas << i << endl;
        auto start1 = std::chrono::high_resolution_clock::now();
        for ( x=0; x<i;x++){
            for(y=0; y<i;y++){
                fWyn[x][y] = (f2dArrayA[x][y]+f2dArrayB[x][y]+f2dArrayC[x][y])*fVector[h][0];
                h++;
                fWynVector[x][0]=fWynVector[x][0] + fWyn[x][y];
            }
            h=0;
        }
        auto end1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time1;

        time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1);
        abcXCzas << time1.count() << endl;

        abcXWyn << "i: "<< i << endl;
        for ( x=0; x<i;x++){
            double out = (double) fWynVector[x][0].get_num() / (double) fWynVector[x][0].get_denum();
            abcXWyn << out << endl;
        }


///////////////////////////A * X//////////////////////
        aXCzas << i << endl;
        auto start2 = std::chrono::high_resolution_clock::now();
        for ( x=0; x<i;x++){
            for(y=0; y<i;y++){
                fWyn2[x][y] = f2dArrayA[x][y]*fVector[h][0];
                h++;
                fWynVector2[x][0]=fWynVector2[x][0]+fWyn2[x][y];
            }
            h=0;
        }
        auto end2 = std::chrono::high_resolution_clock::now();
        double res1Time = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - start2).count();
        //res1Time -= 1000000;
        //res1Time /= 1000000;

        aXCzas << res1Time<< endl;

        aXWyn << "i: "<< i << endl;
        for ( x=0; x<i;x++){
            double out = (double) fWynVector2[x][0].get_num() / (double) fWynVector2[x][0].get_denum();
            aXWyn << out << endl;
        }
///////////////////////////A *(B*C)//////////////////////
        abcCzas << i << endl;
        auto start3 = std::chrono::high_resolution_clock::now();

        for ( x=0; x<i;x++)
            for(y=0; y<i;y++)
                for(m=0;m<i;m++)
                    fWyn3[x][y] = fWyn3[x][y] + (f2dArrayB[x][m] * f2dArrayC[m][y]);

         for ( x=0; x<i;x++)
            for(y=0; y<i;y++)
                for(m=0;m<i;m++){
                    fWynAbc[x][y] = fWynAbc[x][y] + (f2dArrayA[x][m] * fWyn3[m][y]);
                }
        auto end3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time3;

        time3 = std::chrono::duration_cast<std::chrono::duration<double>>(end3 - start3);
        abcCzas << time3.count() << endl;

        abcWyn << "i: "<< i << endl;
        for ( x=0; x<i;x++){
            for(y=0; y<i;y++){
                    double out = (double) fWynAbc[x][y].get_num() / (double) fWynAbc[x][y].get_denum();
                    abcWyn << out;
                    abcWyn << " ";
                }
                abcWyn << endl;
            }
///////////////////////////////////////////////////////////

    }
    abcXCzas.close();
    aXCzas.close();
    abcCzas.close();

    abcXWyn.close();
    aXWyn.close();
    abcWyn.close();
}

template<typename T>
void matrixFD(){
    int i, t, x, y, m, v;  // liczniki
    T num, denum;
    std::ifstream fp("macierze.txt");

    ofstream abcXCzas;
    abcXCzas.open ("FMatrixCzasyABCX.txt");//D
    abcXCzas.precision(PRECD);

    ofstream aXCzas;
    aXCzas.open ("FMatrixCzasyAX.txt");//D
    aXCzas.precision(PRECD);

    ofstream abcCzas;
    abcCzas.open ("FMatrixCzasyABC.txt");//D
    abcCzas.precision(PRECD);

    ofstream abcXWyn;
    abcXWyn.open ("FMatrixWynikiABCX.txt");//D
    abcXWyn.precision(PRECF);//

    ofstream aXWyn;
    aXWyn.open ("FMatrixWynAX.txt");//D
    aXWyn.precision(PRECF);//

    ofstream abcWyn;
    abcWyn.open ("FMatrixWynABC.txt");//D
    abcWyn.precision(PRECF);//

    ofstream fullGaussWyn;
    fullGaussWyn.open ("FMatrixfullGaussWyn.txt");//D
    fullGaussWyn.precision(PRECF);//

    ofstream partGaussWyn;
    partGaussWyn.open ("FMatrixpartGaussWyn.txt");//D
    partGaussWyn.precision(PRECF);//

    ofstream noneGaussWyn;
    noneGaussWyn.open ("FMatrixnoneGaussWyn.txt");//D
    noneGaussWyn.precision(PRECF);//

    ofstream fullGaussCzas;
    fullGaussCzas.open ("FMatrixfullGaussCzas.txt");//D
    fullGaussCzas.precision(PRECD);

    ofstream partGaussCzas;
    partGaussCzas.open ("FMatrixpartGaussCzas.txt");//D
    partGaussCzas.precision(PRECD);

    ofstream noneGaussCzas;
    noneGaussCzas.open ("FMatrixnoneGaussCzas.txt");//D
    noneGaussCzas.precision(PRECD);

    for(i=10;i<MAX;i = i+ 10){
        MyMatrix<T> matrixA(i, i, 0.0);
        MyMatrix<T> matrixB(i, i, 0.0);
        MyMatrix<T> matrixC(i, i, 0.0);
        vector<T> mVector(i);
        vector<T> mVectorWyn(i);
        vector<T> mVectorWyn2(i);
        MyMatrix<T> mWyn(i, i, 0.0);
        MyMatrix<T> mWyn2(i, i, 0.0);
        MyMatrix<T> mWyn3(i, i, 0.0);
        MyMatrix<T> matrixAV(i, i+1, 1.0);
        MyMatrix<T> matrixAV2(i, i+1, 1.0);
        MyMatrix<T> matrixAV3(i, i+1, 1.0);

        for ( x=0; x<i;x++){

            for(y=0; y<i;y++){
                fp >> num;
                fp >> denum;
                T r1 = num / denum;
                matrixA.matrix[x][y] = r1;
                matrixAV.matrix[x][y] = r1;
                matrixAV2.matrix[x][y] = r1;
                matrixAV3.matrix[x][y] = r1;
                fp >> num;
                fp >> denum;
                r1 = num / denum;
                matrixB.matrix[x][y] = r1;
                fp >> num;
                fp >> denum;
                r1 = num / denum;
                matrixC.matrix[x][y] = r1;
            }
            mVector[x] = 0;
            mVectorWyn[x] = 0;
            mVectorWyn2[x] = 0;
        }
        for (t=0; t<i; t++) { // wektor
            T num;
            fp >> num;
            T denum;
            fp >> denum;
            T wyn = num / denum;
            mVector[t] = wyn;
        }
        for(v =0; v<matrixAV.getRowCount(); v++){
                matrixAV(v,matrixAV.getColCount()-1)=mVector[v];
                matrixAV2(v,matrixAV2.getColCount()-1)=mVector[v];
                matrixAV3(v,matrixAV3.getColCount()-1)=mVector[v];
        }

//////////// A * X /////////////
        aXCzas << i << endl;
        auto start1 = std::chrono::high_resolution_clock::now();

        mVectorWyn = matrixA *  mVector;

        auto end1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time1;
        time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1);
        aXCzas << time1.count() << endl;

        aXWyn << "i: "<< i << endl;
        for(x=0; x<i;x++){
            aXWyn << mVectorWyn[x] << endl;
        }
///////////////////////////(A + B + C) * X//////////////////////
        abcXCzas << i << endl;
        auto start2 = std::chrono::high_resolution_clock::now();

        for ( x=0; x<i;x++)
            for(y=0; y<i;y++){
                 mWyn.matrix[x][y] = matrixA.matrix[x][y] + matrixB.matrix[x][y] + matrixC.matrix[x][y];
                }
            mVectorWyn = mWyn * mVector;

        auto end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time2;
        time2 = std::chrono::duration_cast<std::chrono::duration<double>>(end2 - start2);
        abcXCzas << time2.count() << endl;

        abcXWyn << "i: "<< i << endl;
        for(x=0; x<i;x++){
            abcXWyn << mVectorWyn[x] << endl;
        }
///////////////////////////A *(B*C)//////////////////////
        abcCzas << i << endl;
        auto start3 = std::chrono::high_resolution_clock::now();

              for ( x=0; x<i;x++)
                for(y=0; y<i;y++)
                  for(m=0;m<i;m++)
                    mWyn3.matrix[x][y] = mWyn3.matrix[x][y] + (matrixB.matrix[x][m] * matrixC.matrix[m][y]);

              for ( x=0; x<i;x++)
                for(y=0; y<i;y++)
                  for(m=0;m<i;m++)
                    mWyn2.matrix[x][y] = mWyn2.matrix[x][y] + (matrixA.matrix[x][m] * mWyn3.matrix[m][y]);

        auto end3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time3;
        time3 = std::chrono::duration_cast<std::chrono::duration<double>>(end3 - start3);
        abcCzas << time3.count() << endl;

        abcWyn << "i: "<< i << endl;
        for ( x=0; x<i;x++){
            for(y=0; y<i;y++){
                abcWyn << mWyn2.matrix[x][y];
                abcWyn << " ";
            }
            abcWyn << endl;
        }
//////////////////////////////////////////////////////////
        vector<T> gaussFullOut(i);
        vector<T> gaussPartOut(i);
        vector<T> gaussNoneOut(i);
///////////////////////////GAUSS FULL///////////////////////////

        fullGaussCzas << i << endl;
        auto start4 = std::chrono::high_resolution_clock::now();

        gaussFullOut = matrixAV.Gauss1();

        auto end4 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time4;
        time4 = std::chrono::duration_cast<std::chrono::duration<double>>(end4 - start4);
        fullGaussCzas << time4.count() << endl;

        fullGaussWyn << i << endl;
        for ( x=0; x<i;x++)
            fullGaussWyn << gaussFullOut[x] << endl;

/////////////////////////GAUSS PARTIAL/////////////////////////
        partGaussCzas << i << endl;
        auto start5 = std::chrono::high_resolution_clock::now();

        gaussPartOut = matrixAV2.Gauss2();

        auto end5 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time5;
        time5 = std::chrono::duration_cast<std::chrono::duration<double>>(end5 - start5);
        partGaussCzas << time5.count() << endl;

        partGaussWyn << i << endl;
        for ( x=0; x<i;x++)
            partGaussWyn << gaussPartOut[x] << endl;


////////////////////////////GAUSS NONE/////////////////////////////
        noneGaussCzas << i << endl;
        auto start6 = std::chrono::high_resolution_clock::now();

        gaussNoneOut = matrixAV3.Gauss3();

        auto end6 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time6;
        time6 = std::chrono::duration_cast<std::chrono::duration<double>>(end6 - start6);
        noneGaussCzas << time6.count() << endl;

        noneGaussWyn << i << endl;
        for ( x=0; x<i;x++)
            noneGaussWyn << gaussNoneOut[x] << endl;

        }
    abcXCzas.close();
    aXCzas.close();
    abcCzas.close();

    abcXWyn.close();
    aXWyn.close();
    abcWyn.close();

    fullGaussWyn.close();
    noneGaussWyn.close();
    partGaussWyn.close();

    fullGaussCzas.close();
    noneGaussCzas.close();
    partGaussCzas.close();
};

void eigenFloat(){
    int i, t, x, y;
    float num, denum;
    std::ifstream fp("macierze.txt");
    ofstream abcXCzas;
    abcXCzas.open ("FEigenCzasyABCX.txt");
    abcXCzas.precision(PRECD);

    ofstream aXCzas;
    aXCzas.open ("FEigenCzasyAX.txt");
    aXCzas.precision(PRECD);

    ofstream abcCzas;
    abcCzas.open ("FEigenCzasyABC.txt");
    abcCzas.precision(PRECD);

    ofstream abcXWyn;
    abcXWyn.open ("FEigenWynikiABCX.txt");
    abcXWyn.precision(PRECF);

    ofstream aXWyn;
    aXWyn.open ("FEigenWynAX.txt");
    aXWyn.precision(PRECF);

    ofstream abcWyn;
    abcWyn.open ("FEigenWynABC.txt");
    abcWyn.precision(PRECF);

    ofstream fullGaussWyn;
    fullGaussWyn.open ("FEigenfullGaussWyn.txt");
    fullGaussWyn.precision(PRECF);

    ofstream partGaussWyn;
    partGaussWyn.open ("FEigenpartGaussWyn.txt");
    partGaussWyn.precision(PRECF);

    ofstream fullGaussCzas;
    fullGaussCzas.open ("FEigenfullGaussCzas.txt");
    fullGaussCzas.precision(PRECD);

    ofstream partGaussCzas;
    partGaussCzas.open ("FEigenpartGaussCzas.txt");
    partGaussCzas.precision(PRECD);


        for(i=10;i<MAX;i = i+ 10){
        MatrixXf emA(i,i);
        MatrixXf emB(i,i);
        MatrixXf emC(i,i);
        VectorXf evV(i);
        VectorXf evW(i);
        MatrixXf emW(i,i);
        MatrixXf emW2(i,i);

        for ( x=0; x<i;x++){
            evW(x) = 0.0;
            for(y=0; y<i;y++){
                emW(x,y) = 0.0;
                emW2(x,y) = 0.0;
                fp >> num;
                fp >> denum;
                float r1 = num / denum;
                emA(x,y) = r1;
                fp >> num;
                fp >> denum;
                r1 = num / denum;
                emB(x,y) = r1;
                fp >> num;
                fp >> denum;
                r1 = num / denum;
                emC(x,y) = r1;
            }
        }

        for (t=0; t<i; t++) { // wektor
            float num;
            fp >> num;
            float denum;
            fp >> denum;
            float wyn = num / denum;
            evV(t) = wyn;
        }

///////////////// A * X /////////////////
        aXCzas << i << endl;
        auto start1 = std::chrono::high_resolution_clock::now();

        evW = emA * evV;

        auto end1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time1;

        time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1);
        aXCzas << time1.count() << endl;

        aXWyn << "i: "<< i << endl;
        aXWyn << evW << endl;
///////////////////////////(A + B + C) * X//////////////////////
        abcXCzas << i << endl;
        auto start2 = std::chrono::high_resolution_clock::now();

        evW = (emA + emB + emC) * evV;

        auto end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time2;

        time2 = std::chrono::duration_cast<std::chrono::duration<double>>(end2 - start2);
        abcXCzas << time2.count() << endl;

        abcXWyn << "i: "<< i << endl;
        abcXWyn << evW << endl;
///////////////////////////A *(B*C)//////////////////////
        abcCzas << i << endl;
        auto start3 = std::chrono::high_resolution_clock::now();

        emW = emB * emC;
        emW2 = emA * emW;
        //std::cout << emB*emC;
        auto end3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time3;

        time3 = std::chrono::duration_cast<std::chrono::duration<double>>(end3 - start3);
        abcCzas << time3.count() << endl;

        abcWyn << "i: "<< i << endl;
        abcWyn << emW2 << endl;
//////////////////////////////////////////////////////////
                    // fullPivLu
        fullGaussCzas << i << endl;
        auto start4 = std::chrono::high_resolution_clock::now();

        VectorXf fullGauss = emA.fullPivLu().solve(evV);

        auto end4 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time4;
        time4 = std::chrono::duration_cast<std::chrono::duration<double>>(end4 - start4);
        fullGaussCzas << time4.count() << endl;

        fullGaussWyn << i << endl;
        fullGaussWyn << fullGauss << endl;

            // PartialPivLu
        auto start5 = std::chrono::high_resolution_clock::now();

        VectorXf partialGauss = emA.partialPivLu().solve(evV);

        auto end5 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time5;
        time5 = std::chrono::duration_cast<std::chrono::duration<double>>(end5 - start5);
        partGaussCzas << time5.count() << endl;

        partGaussWyn << i << endl;
        partGaussWyn << partialGauss << endl;
    }
    abcXCzas.close();
    aXCzas.close();
    abcCzas.close();

    abcXWyn.close();
    aXWyn.close();
    abcWyn.close();

    fullGaussWyn.close();
    partGaussWyn.close();

    fullGaussCzas.close();
    partGaussCzas.close();
}

void eigenDouble(){
    int i, t, x, y;
    double num, denum;
    std::ifstream fp("macierze.txt");
    ofstream abcXCzas;
    abcXCzas.open ("DEigenCzasyABCX.txt");
    abcXCzas.precision(PRECD);

    ofstream aXCzas;
    aXCzas.open ("DEigenCzasyAX.txt");
    aXCzas.precision(PRECD);

    ofstream abcCzas;
    abcCzas.open ("DEigenCzasyABC.txt");
    abcCzas.precision(PRECD);

    ofstream abcXWyn;
    abcXWyn.open ("DEigenWynikiABCX.txt");
    abcXWyn.precision(PRECD);

    ofstream aXWyn;
    aXWyn.open ("DEigenWynAX.txt");
    aXWyn.precision(PRECD);

    ofstream abcWyn;
    abcWyn.open ("DEigenWynABC.txt");
    abcWyn.precision(PRECD);

    ofstream fullGaussWyn;
    fullGaussWyn.open ("DEigenfullGaussWyn.txt");
    fullGaussWyn.precision(PRECD);

    ofstream partGaussWyn;
    partGaussWyn.open ("DEigenpartGaussWyn.txt");
    partGaussWyn.precision(PRECD);

    ofstream fullGaussCzas;
    fullGaussCzas.open ("DEigenfullGaussCzas.txt");
    fullGaussCzas.precision(PRECD);

    ofstream partGaussCzas;
    partGaussCzas.open ("DEigenpartGaussCzas.txt");
    partGaussCzas.precision(PRECD);

        for(i=10;i<MAX;i = i+ 10){
        MatrixXd emA(i,i);
        MatrixXd emB(i,i);
        MatrixXd emC(i,i);
        VectorXd evV(i);
        VectorXd evW(i);
        MatrixXd emW(i,i);
        MatrixXd emW2(i,i);

        for ( x=0; x<i;x++){
            evW(x) = 0.0;
            for(y=0; y<i;y++){
                emW(x,y) = 0.0;
                emW2(x,y) = 0.0;
                fp >> num;
                fp >> denum;
                double r1 = num / denum;
                emA(x,y) = r1;
                fp >> num;
                fp >> denum;
                r1 = num / denum;
                emB(x,y) = r1;
                fp >> num;
                fp >> denum;
                r1 = num / denum;
                emC(x,y) = r1;
            }
        }

        for (t=0; t<i; t++) { // wektor
            double num;
            fp >> num;
            double denum;
            fp >> denum;
            double wyn = num / denum;
            evV(t) = wyn;
        }

///////////////// A * X /////////////////
        aXCzas << i << endl;
        auto start1 = std::chrono::high_resolution_clock::now();

        evW = emA * evV;

        auto end1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time1;

        time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1);
        aXCzas << time1.count() << endl;

        aXWyn << "i: "<< i << endl;
        aXWyn << evW << endl;
///////////////////////////(A + B + C) * X//////////////////////
        abcXCzas << i << endl;
        auto start2 = std::chrono::high_resolution_clock::now();

        evW = (emA + emB + emC) * evV;

        auto end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time2;

        time2 = std::chrono::duration_cast<std::chrono::duration<double>>(end2 - start2);
        abcXCzas << time2.count() << endl;

        abcXWyn << "i: "<< i << endl;
        abcXWyn << evW << endl;
///////////////////////////A *(B*C)//////////////////////
        abcCzas << i << endl;
        auto start3 = std::chrono::high_resolution_clock::now();

        emW = emB * emC;
        emW2 = emA * emW;

        auto end3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time3;

        time3 = std::chrono::duration_cast<std::chrono::duration<double>>(end3 - start3);
        abcCzas << time3.count() << endl;

        abcWyn << "i: "<< i << endl;
        abcWyn << emW2 << endl;
/////////////////////////////////////////////////////////

                    // fullPivLu
        fullGaussCzas << i << endl;
        auto start4 = std::chrono::high_resolution_clock::now();

        VectorXd fullGauss = emA.fullPivLu().solve(evV);

        auto end4 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time4;
        time4 = std::chrono::duration_cast<std::chrono::duration<double>>(end4 - start4);
        fullGaussCzas << time4.count() << endl;

        fullGaussWyn << i << endl;
        fullGaussWyn << fullGauss << endl;

            // PartialPivLu
        auto start5 = std::chrono::high_resolution_clock::now();

        VectorXd partialGauss = emA.partialPivLu().solve(evV);

        auto end5 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time5;
        time5 = std::chrono::duration_cast<std::chrono::duration<double>>(end5 - start5);
        partGaussCzas << time5.count() << endl;

        partGaussWyn << i << endl;
        partGaussWyn << partialGauss << endl;

    }
    abcXCzas.close();
    aXCzas.close();
    abcCzas.close();

    abcXWyn.close();
    aXWyn.close();
    abcWyn.close();

    fullGaussWyn.close();
    partGaussWyn.close();

    fullGaussCzas.close();
    partGaussCzas.close();
}

void matrixNorm(){
    int x,t,y;
    std::ifstream fMf("FMatrixfullGaussWyn.txt");
    std::ifstream fMp("FMatrixpartGaussWyn.txt");
    std::ifstream fMn("FMatrixnoneGaussWyn.txt");
  /*std::ifstream dMf("DMatrixfullGaussWyn.txt");
    std::ifstream dMp("DMatrixpartGaussWyn.txt");
    std::ifstream dMn("DMatrixnoneGaussWyn.txt");*/

    std::ifstream fEf("FEigenfullGaussWyn.txt");
    std::ifstream fEp("FEigenpartGaussWyn.txt");
  //std::ifstream dEf("DEigenfullGaussWyn.txt");
  //std::ifstream dEp("DEigenpartGaussWyn.txt");


    ofstream FfullGaussNorm;
    FfullGaussNorm.open ("FfullGaussNorm.txt");
    FfullGaussNorm.precision(PRECF);
    ofstream FpartGaussNorm;
    FpartGaussNorm.open ("FpartGaussNorm.txt");
    FpartGaussNorm.precision(PRECF);
    ofstream FnoneGaussNorm;
    FnoneGaussNorm.open ("FnoneGaussNorm.txt");
    FnoneGaussNorm.precision(PRECF);
  /*ofstream DfullGaussNorm;
    DfullGaussNorm.open ("DfullGaussNorm.txt");
    DfullGaussNorm.precision(PRECF);
    ofstream DpartGaussNorm;
    DpartGaussNorm.open ("DpartGaussNorm.txt");
    DpartGaussNorm.precision(PRECF);
    ofstream DnoneGaussNorm;
    DnoneGaussNorm.open ("DnoneGaussNorm.txt");
    DnoneGaussNorm.precision(PRECF);*/


    for(x=10;x<MAX;x=x+10){
        float ef, efF, mf, xff=0.0, xfp=0.0, xfn=0.0;
        double ed, edF, md, xdf=0.0, xdp=0.0, xdn=0.0;
        for(y=0;y<x;y++){
            fMf >> t; fMp >> t; fMn >> t;
          //dMf >> t; dMp >> t; dMn >> t; 

            fEf >> t; fEp >> t;
          //dEf >> t; dEp >> t;

            fMf >> mf;
            fEf >> efF;
            xff += pow((mf-efF), 2);

            fMp >> mf;
            fEp >> ef;
            xfp += pow((mf-ef), 2);

            fMn >> mf;
            xfn += pow((mf-efF), 2);
          /*dMf >> md;
            dEf >> edf;
            xdf += pow((md-edf),2);

            dMp >> md;
            dEp >> ed;
            xdp += pow((md-ed),2);

            dMn >> md;
            xdn += pow((md-edf), 2);*/
        }
        FfullGaussNorm << sqrt(xff) << endl;
        FpartGaussNorm << sqrt(xfp) << endl;
        FnoneGaussNorm << sqrt(xfn) << endl;
      /*DfullGaussNorm << sqrt(xdf) << endl;
        DpartGaussNorm << sqrt(xdp) << endl;
        DnoneGaussNorm << sqrt(xdn) << endl;*/
    }
    FfullGaussNorm.close();
    FpartGaussNorm.close();
    FnoneGaussNorm.close();
  /*DfullGaussNorm.close();
    DpartGaussNorm.close();
    DnoneGaussNorm.close();*/
}

int main() {
    //matrixFD<double>();
    matrixFD<float>();
    cout << "matrix float exit \n";
    eigenDouble();
    cout << "eigen double exit \n";
    eigenFloat();
    cout << "eigen float exit \n";
    ulamek();
    cout << "ulamek exit \n";

    matrixNorm();

    int quit;
    scanf("%d", &quit);
    return 0;
}
