#ifndef MATLIB_H
#define MATLIB_H
#include <string>
#include <cmath>

namespace Matlib{ //Declaration of classes
    template<int r, int c>
    class Matrix;

    class vec3;
    class vec4;
    template<int n>
    class vec;
}

namespace Matlib{
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

        Matrix add(Matrix const &m);
        Matrix operator+(Matrix const &m);

        template<int t>
        Matrix<r,t> mult(const Matrix<c,t> &m);
        template<int t>
        Matrix<r,t> operator* (const Matrix<c,t> &m);
        
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
        std::string size();

        Matrix copy();

        Matrix<n,n> transpose();

        Matrix add(Matrix const &m);
        Matrix operator+(Matrix const &m);

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
        Matrix getInverse();

        private:
        int getTriangle(Matrix* matOut);
    };
}

namespace Matlib{
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

#endif //MATLIB_H Header Guard