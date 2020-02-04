#include "matlib.hpp"
#define FLOAT_TOLERANCE 0.000000001

namespace Matlib{ //NxM Matrix Definitions
    template<int r, int c>
    Matrix<r,c>::Matrix(){ //Default Constructor
        double vals[r][c] = {0};
        memcpy(elems, vals, sizeof(double)*r*c);
    }

    template<int r, int c>
    Matrix<r,c>::Matrix(double vals[r][c]){ //Main constructor from 2D array
        memcpy(elems, vals, sizeof(double)*r*c);
    }

    template<int r, int c>
    Matrix<r,c>::Matrix(const Matrix &m){
        memcpy(elems, m.elems, sizeof(double)*r*c);
    }

    template<int r, int c>
    Matrix<r,c> Matrix<r,c>::operator= (const Matrix &m){
            memcpy(elems, m.elems, sizeof(double)*r*c);
    }

    template<int r, int c>
    double* Matrix<r,c>::operator [](int i){
        if(i < r){
            return elems[i];
        }
        return nullptr;
    }

    template<int r, int c>
    std::string Matrix<r,c>::to_string(){
        std::string s = "";
        for(int i = 0; i < r; i++){
            for(int j = 0; j < c; j++){
                s += std::to_string(elems[i][j]);
                if(j < c-1) s+=" ";
            }
            if(i < r-1)
                s += '\n';
        }
        return s;
    }

    template<int r, int c>
    typename Matrix<r,c>::Size Matrix<r,c>::size(){
        return Size(r,c);
    }

    template<int r, int c>
    Matrix<r,c> Matrix<r,c>::copy(){
        return Matrix<r,c>(*this);
    }

    template<int r, int c>
    Matrix<c,r> Matrix<r,c>::transpose(){
        double vals[c][r];
        for(int i = 0; i < c; i++){
            for(int j = 0; j < r; j++){
                vals[i][j] = elems[j][i];
            }
        }
        return Matrix<c,r>(vals);
    }

    template<int r, int c>
    Matrix<r,c> Matrix<r,c>::add(Matrix const &m){
        Matrix<r,c> result = copy();
        for(int i = 0; i < r; i++){
            for(int j = 0; j < c; j++){
                result[i][j] += m[i][j];
            }
        }
        return result;
    }

    template<int r, int c>
    Matrix<r,c> Matrix<r,c>::operator+(Matrix const &m){
        return add(m);
    }

    template<int r, int c>
    template<int t>
    Matrix<r,t> Matrix<r,c>::mult(const Matrix<c,t> &m){
        double vals[r][t];
        for(int i = 0; i < r; i++){
            for(int j = 0; j < t; j++){
                double val = 0;
                for(int k = 0; k < c; k++){
                    val += elems[i][k]*m[k][j];
                }
                
                vals[i][j] = val;
            }
        }
        Matrix<r,t> result(vals);
        return result;
    }

    template<int r, int c>
    template<int t>
    Matrix<r,t> Matrix<r,c>::operator*(const Matrix<c,t> &m){
        return mult(m);
    }

    template<int r, int c>
    Matrix<r,c> Matrix<r,c>::mult(double s){
        Matrix m = copy();
        for(int i = 0; i < r; i++){
            for(int j = 0; j < c; j++){
                m[i][j] *= s;
            }
        }
        return m;
    }

    template<int r, int c>
    Matrix<r,c> Matrix<r,c>::operator*(double s){
        return mult(s);
    }

}

namespace Matlib{//NxN Matrix Definitions
    template<int n>
    Matrix<n,n>::Matrix(){
        double vals[n][n] = {0};
        memcpy(elems, vals, sizeof(double)*n*n);
    }
}