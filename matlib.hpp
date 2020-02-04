#ifndef MATLIB_H
#define MATLIB_H
#include <string>
#include <cmath>

namespace matlib{ //Definition of classes
    template<int r, int c>
    class matrix;

    class vec3;
    class vec4;
    template<int n>
    class vec;
}

namespace matlib{
    template<int r, int c>
    class matrix{
        public:
        double elems[r][c];

        matrix(double vals[r][c]);
        matrix(const matrix &m);

        matrix operator= (const matrix &m); //For assignment to other variables
        matrix operator= (const matrix m); //For assignment to evaluated matrices

        double* operator [](int i);

        std::string to_string();

        std::string size();

        matrix copy();

        matrix<c,r> transpose();

        matrix add(matrix const &m);

        matrix operator+(matrix const &m);

        template<int t>
        matrix<r,t> mult(matrix<c,t> &m);
        template<int t>
        matrix<r,t> operator* (matrix<c,t> &m);
        
        matrix mult(double s);
        matrix operator* (double s);

        private:
        matrix getLowerTriangle();
        matrix getUpperTriangle();
    };

    template<int n>
    class matrix<n,n>{
        public:
        double elems[n][n];

        matrix(double vals[n][n]);

        matrix(const matrix &m);

        matrix operator= (const matrix &m);

        matrix operator= (const matrix m);

        static matrix Identity();

        double* operator [](int i);

        std::string to_string();

        std::string size();

        matrix copy();

        matrix<n,n> transpose();

        matrix add(matrix const &m);

        matrix operator+(matrix const &m);

        template<int t>
        matrix<n,t> mult(matrix<n,t> &m);
        template<int t>
        matrix<n,t> operator* (matrix<n,t> &m);

        matrix mult(double s);
        matrix operator* (double s);

        double determinant();

        matrix getMinors();
        matrix getCofactor();
        matrix getAdjugate();
        matrix getInverse();
    };
}

namespace matlib{
    class vec3{
        public:
        double x,y,z;
        
        vec3();
        vec3(double x, double y, double z);
        vec3(matrix<3,1> m);

        matrix<3,1> getMatrix();
        std::string to_string();

        double mag();

        double dot(vec3 v);
        vec3 cross(vec3 v);

        double angleBetween(vec3 v);


        vec3 transform(matrix<3,3> &m);
    };
    class vec4{
        public:
        double x,y,z,w;
        vec4(double x, double y, double z);
        vec4(vec3 v);
        vec4(matrix<4,1> m);

        matrix<4,1> getMatrix();

        vec4 transform(matrix<4,4> &m);
    };
}

#endif //MATLIB_H Header Guard