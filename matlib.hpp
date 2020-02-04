#ifndef MATLIB_H
#define MATLIB_H
#include <string>
#include <cmath>

namespace matlib{
    template<int r, int c, class ElemType=double>
    class matrix{
        public:
        ElemType elems[r][c];

        matrix(ElemType vals[r][c]);
        matrix(const matrix &m);

        matrix operator= (const matrix &m); //For assignment to other variables
        matrix operator= (const matrix m); //For assignment to evaluated matrices

        ElemType* operator [](int i);

        std::string to_string();

        std::string size();

        matrix copy();

        matrix<c,r> transpose();

        matrix add(matrix const &m);

        matrix operator+(matrix const &m);

        template<int t>
        matrix<r,t,ElemType> mult(matrix<c,t,ElemType> &m);
        template<int t>
        matrix<r,t,ElemType> operator* (matrix<c,t,ElemType> &m);
        
        matrix mult(ElemType s);
        matrix operator* (ElemType s);

        private:
        matrix getLowerTriangle();
        matrix getUpperTriangle();
    };

    template<int n, class ElemType>
    class matrix<n,n,ElemType>{
        public:
        ElemType elems[n][n];

        matrix(ElemType vals[n][n]);

        matrix(const matrix &m);

        matrix operator= (const matrix &m);

        matrix operator= (const matrix m);

        static matrix Identity();

        ElemType* operator [](int i);

        std::string to_string();

        std::string size();

        matrix copy();

        matrix<n,n> transpose();

        matrix add(matrix const &m);

        matrix operator+(matrix const &m);

        template<int t>
        matrix<n,t,ElemType> mult(matrix<n,t,ElemType> &m);
        template<int t>
        matrix<n,t> operator* (matrix<n,t,ElemType> &m);

        matrix mult(ElemType s);
        matrix operator* (ElemType s);

        ElemType determinant();

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