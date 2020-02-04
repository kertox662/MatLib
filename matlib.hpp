#ifndef MATLIB_H //MATLIB_H Header Guard
#define MATLIB_H

#include <string>
#include <cmath>
#define MATLIB_FLOAT_TOLERANCE 0.000000001

namespace Matlib{ //Declaration of classes
    template<int r, int c>
    class Matrix;

    class vec3;
    class vec4;
    template<int n>
    class vec;
}

namespace Matlib{ //Matrix Class Member Declerations
    template<int r, int c>
    class Matrix{
        double elems[r][c];

        public:
        Matrix();
        Matrix(double vals[r][c]);

        Matrix(const Matrix &m);
        Matrix operator= (const Matrix &m);

        double* operator [](int i);

        std::string to_string();

        struct Size{
            int rows,columns;
            std::string to_string(){
                return std::to_string(rows) + 'x' + std::to_string(columns);
            }
        };
        
        Size size();

        Matrix copy();

        Matrix<c,r> transpose();

        Matrix add(Matrix &m);
        Matrix operator+(Matrix &m);

        template<int t>
        Matrix<r,t> mult(Matrix<c,t> &m);
        template<int t>
        Matrix<r,t> operator* (Matrix<c,t> &m);
        
        Matrix mult(double s);
        Matrix operator* (double s);

    };

    template<int n>
    class Matrix<n,n>{
        double elems[n][n];

        public:
        Matrix();
        Matrix(double vals[n][n]);

        Matrix(const Matrix &m);
        Matrix operator= (const Matrix &m);

        static Matrix Identity();

        double* operator [](int i);

        std::string to_string();

        struct Size{
            int rows,columns;
            std::string to_string(){
                return std::to_string(rows) + 'x' + std::to_string(columns);
            }
        };

        Size size();

        Matrix copy();

        Matrix<n,n> transpose();

        Matrix add(Matrix &m);
        Matrix operator+(Matrix &m);

        template<int t>
        Matrix<n,t> mult(Matrix<n,t> &m);
        template<int t>
        Matrix<n,t> operator* (Matrix<n,t> &m);

        Matrix mult(double s);
        Matrix operator* (double s);

        double determinant();

        Matrix getMinors();
        Matrix getCofactor();
        Matrix getAdjugate();
        Matrix Inverse();
        Matrix Inverse(int* error);

        private:
        int getTriangles(Matrix* upperOut){
            Matrix upper = copy();
            int P[n]; //Permutation matrix;
            int numPivots = 0;
            for(int i = 0; i < n; i++)
                P[i]=i;
            for(int i = 0; i < n; i++){

                double maxV = 0, absV;
                int maxInd = i;

                for(int k = i; k < n; k++){
                    if((absV = fabs(upper[k][i])) > maxV){
                        maxV = absV;
                        maxInd = k;
                    }
                }

                if (maxV < MATLIB_FLOAT_TOLERANCE) return -1; //Degenerate Matrix, Determinant 0

                if(maxInd != i){
                    int permSwap = P[i]; //Permuation matrix swapper
                    P[i] = P[maxInd];
                    P[maxInd] = P[permSwap];

                    for(int c = 0; c < n; c++){
                        double rowSwap = upper[i][c];
                        upper[i][c] = upper[maxInd][c];
                        upper[maxInd][c] = rowSwap;
                    }

                    numPivots++;
                }

                for(int j = i+1; j < n; j++){
                    double factor = upper[j][i] / upper[i][i];
                    for(int k = 0; k < n; k++){
                        upper[j][k] -= factor*upper[i][k];
                        
                    }
                }
            }
            *upperOut = upper;

            return numPivots;
        }
    };
}

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
        return copy();
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
    Matrix<r,c> Matrix<r,c>::add(Matrix &m){
        Matrix result = copy();
        for(int i = 0; i < r; i++){
            for(int j = 0; j < c; j++){
                result[i][j] += m[i][j];
            }
        }
        return result;
    }

    template<int r, int c>
    Matrix<r,c> Matrix<r,c>::operator+(Matrix &m){
        return add(m);
    }

    template<int r, int c>
    template<int t>
    Matrix<r,t> Matrix<r,c>::mult(Matrix<c,t> &m){
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
    Matrix<r,t> Matrix<r,c>::operator*(Matrix<c,t> &m){
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

namespace Matlib{ //NxN Matrix Definitions
    template<int n>
    Matrix<n,n>::Matrix(){ //Default Constructor
        double vals[n][n] = {0};
        memcpy(elems, vals, sizeof(double)*n*n);
    }

    template<int n>
    Matrix<n,n>::Matrix(double vals[n][n]){ //Main constructor from 2D array
        memcpy(elems, vals, sizeof(double)*n*n);
    }

    template<int n>
    Matrix<n,n>::Matrix(const Matrix &m){
        memcpy(elems, m.elems, sizeof(double)*n*n);
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::operator= (const Matrix &m){
        memcpy(elems, m.elems, sizeof(double)*n*n);
        return copy();
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::Identity(){
        Matrix m;
        for(int i = 0; i < n; i++){
            m[i][i] = 1;
        }
        return m;
    }

    template<int n>
    double* Matrix<n,n>::operator [](int i){
        if(i < n){
            return elems[i];
        }
        return nullptr;
    }

    template<int n>
    std::string Matrix<n,n>::to_string(){
        std::string s = "";
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                s += std::to_string(elems[i][j]);
                if(j < n-1) s+=" ";
            }
            if(i < n-1)
                s += '\n';
        }
        return s;
    }

    template<int n>
    typename Matrix<n,n>::Size Matrix<n,n>::size(){
        return Size(n,n);
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::copy(){
        return Matrix<n,n>(*this);
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::transpose(){
        double vals[n][n];
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                vals[i][j] = elems[j][i];
            }
        }
        return Matrix<n,n>(vals);
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::add(Matrix &m){
        Matrix result = copy();
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                result[i][j] += m[i][j];
            }
        }
        return result;
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::operator+(Matrix &m){
        return add(m);
    }

    template<int n>
    template<int t>
    Matrix<n,t> Matrix<n,n>::mult(Matrix<n,t> &m){
        double vals[n][t];
        for(int i = 0; i < n; i++){
            for(int j = 0; j < t; j++){
                double val = 0;
                for(int k = 0; k < n; k++){
                    val += elems[i][k]*m[k][j];
                }
                
                vals[i][j] = val;
            }
        }
        return Matrix<n,t>(vals);
    }

    template<int n>
    template<int t>
    Matrix<n,t> Matrix<n,n>::operator*(Matrix<n,t> &m){
        return mult(m);
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::mult(double s){
        Matrix m = copy();
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                m[i][j] *= s;
            }
        }
        return m;
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::operator*(double s){
        return mult(s);
    }

    template<int n>
    double Matrix<n,n>::determinant(){
        Matrix m;
        int nPivot = getTriangles(&m);
        if(nPivot == -1) return 0;
        int d = 1;
        for(int i = 0; i < n; i++){
            d *= m[i][i];
        }
        if(nPivot%2)
            return -d;
        else
            return d;
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::getMinors(){
        Matrix minors;
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                double vals[n-1][n-1];
                int indR = 0, indC = 0;
                for(int r = 0; r < n; r++){
                    for(int c = 0; c < n; c++){
                        if(r != i && c != j){
                            vals[indR][indC] = elems[r][c];
                            indC++;
                            if(indC >= n-1){
                                indC = 0;
                                indR++;
                            }
                        }
                    }
                }
                Matrix<n-1,n-1> minor(vals);
                minors[i][j] = minor.determinant();
            }
        }
        return minors;
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::getCofactor(){
        Matrix minors = getMinors();
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if((i+j)%2){
                    minors[i][j]*=-1;
                }
            }
        }
        return minors;
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::getAdjugate(){
        return getCofactor().transpose();
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::Inverse(){
        double det = determinant();
        if(fabs(det) < MATLIB_FLOAT_TOLERANCE){
            return Matrix();
        }
        return getAdjugate() * (1/det);
    }

    template<int n>
    Matrix<n,n> Matrix<n,n>::Inverse(int* error){
        double det = determinant();
        if(fabs(det) < MATLIB_FLOAT_TOLERANCE){
            *error = 1;
            return Matrix();
        }
        *error = 0;
        return getAdjugate() * (1/det);
    }

}

namespace Matlib{ //Vector Class Member Declaration and Definitions
    class vec3{
        public:
        double x,y,z;
        
        vec3();
        vec3(double x, double y, double z);
        vec3(Matrix<3,1> m);

        Matrix<3,1> getMatrix();
        std::string to_string();

        double mag();

        double dot(vec3 v);
        vec3 cross(vec3 v);

        double angleBetween(vec3 v);


        vec3 transform(Matrix<3,3> &m);
    };
    class vec4{
        public:
        double x,y,z,w;
        vec4(double x, double y, double z);
        vec4(vec3 v);
        vec4(Matrix<4,1> m);

        Matrix<4,1> getMatrix();

        vec4 transform(Matrix<4,4> &m);
    };
}

namespace Matlib{ //Functions for various conversions
    namespace Conversion{
        double rad_to_deg(double r) {
            return r * (180.0 / M_PI);
        }

        double deg_to_rad(double d){
            return d * M_PI / 180.0;
        }
    }
}

namespace Matlib{ //Functions that help with transformations
    namespace Transform{
        Matrix<2,2> make2DRot(double a){
            double v[2][2] = {{cos(a),-sin(a)}, {sin(a), cos(a)}};
            return Matrix<2,2>(v);
        }

        Matrix<3,3> makeXRot(double a){
            double vals[3][3] = {{1,0,0},{0,cos(a),-sin(a)},{0,sin(a),cos(a)}};
            return Matrix<3,3>(vals);
        }

        Matrix<3,3> makeYRot(double a){
            double vals[3][3] = {{cos(a),0,sin(a)},{0,1,0},{-sin(a),0,cos(a)}};
            return Matrix<3,3>(vals);
        }

        Matrix<3,3> makeZRot(double a){
            double vals[3][3] = {{cos(a),-sin(a),0},{sin(a),cos(a),0},{0,0,1}};
            return Matrix<3,3>(vals);
        }

        Matrix<3,3> makeRot(double a, double b, double c){
            Matrix<3,3> mx = makeXRot(a), my = makeYRot(b), mz = makeZRot(c);
            return mx * my * mz;
        }

        Matrix<4,4> makeTranslation(vec3 v){
            double vals[4][4] = {{1,0,0,v.x},{0,1,0,v.y},{0,0,1,v.z},{0,0,0,1}};
            return Matrix<4,4>(vals);
        }

        Matrix<3,3> makeScale(double scale){
            return Matrix<3,3>::Identity() * scale;
        }
    }
}

#endif //MATLIB_H Header Guard