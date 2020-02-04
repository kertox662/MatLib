// #include "matlib.hpp"
#include <string>
#include <cmath>
#include <iostream>
#define FLOAT_TOLERANCE 0.000000001

using std::cout;
using std::endl;

namespace matlib{
    template<int r, int c>
    class matrix{
        public:
        double elems[r][c];

        matrix(){
            double vals[r][c] = {0};
            memcpy(elems, vals, sizeof(double)*r*c);
        }

        matrix(double vals[r][c]){
            memcpy(elems, vals, sizeof(double)*r*c);
        }

        matrix(const matrix &m){
            memcpy(elems, m.elems, sizeof(double)*r*c);
        }

        matrix operator= (const matrix &m){
            memcpy(elems, m.elems, sizeof(double)*r*c);
        }

        double* operator [](int i){
            if(i < r){
                return elems[i];
            }
            return nullptr;
        }

        std::string to_string(){
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

        std::string size(){
            return to_string(r) + 'x' + to_string(c);
        }

        matrix copy(){
            matrix m(*this);
            return m;
        }

        matrix<c,r> transpose(){
            double vals[c][r];
            for(int i = 0; i < c; i++){
                for(int j = 0; j < r; j++){
                    vals[i][j] = elems[j][i];
                }
            }
            matrix<c,r> m(vals);
            return m;
        }

        matrix add(matrix const &m){
            matrix result(elems);
            for(int i = 0; i < r; i++){
                for(int j = 0; j < c; j++){
                    result[i][j] += m[i][j];
                }
            }
            return result;
        }

        matrix operator+(matrix const &m){
            return add(m);
        }

        template<int t>
        matrix<r,t> mult(matrix<c,t> &m){
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
            matrix<r,t> result(vals);
            return result;
        }

        matrix mult(double s){
            matrix m = copy();
            for(int i = 0; i < r; i++){
                for(int j = 0; j < c; j++){
                    m[i][j] *= s;
                }
            }
            return m;
        }

        template<int t>
        matrix<r,t> operator* (matrix<c,t> &m){
            return mult(m);
        }

        matrix operator* (double s){
            return mult(s);
        }
    };

    template<int n>
    class matrix<n,n>{
        public:
        double elems[n][n];

        matrix(){
            double vals[n][n] = {0};
            memcpy(elems, vals, sizeof(double)*n*n);
        }

        matrix(double vals[n][n]){
            memcpy(elems, vals, sizeof(double)*n*n);
        }

        matrix(const matrix &m){
            memcpy(elems, m.elems, sizeof(double) * n * n);
        }

        void operator= (const matrix &m){
            memcpy(elems, m.elems, sizeof(double)* n * n);
        }

        static matrix Identity(){
            double vals[n][n] = {0};
            for(int i = 0; i < n; i++){
                vals[i][i] = 1;
            }
            return matrix(vals);
        }

        double* operator [](int i){
            if(i < n){
                return elems[i];
            }
            return nullptr;
        }

        std::string to_string(){
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

        std::string size(){
            return to_string(n) + 'x' + to_string(n);
        }

        matrix copy(){
            matrix m(*this);
            return m;
        }

        matrix<n,n> transpose(){
            double vals[n][n];
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    vals[i][j] = elems[j][i];
                }
            }
            matrix m(vals);
            return m;
        }

        matrix add(matrix const &m){
            matrix result(elems);
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    result[i][j] += m[i][j];
                }
            }
            return result;
        }

        matrix operator+(matrix const &m){
            return add(m);
        }

        template<int t>
        matrix<n,t> mult(matrix<n,t> &m){
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
            matrix<n,t> result(vals);
            return result;
        }

        matrix mult(double s){
            matrix m = copy();
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    m[i][j] *= s;
                }
            }
            return m;
        }

        template<int t>
        matrix<n,t> operator* (matrix<n,t> &m){
            return mult(m);
        }

        matrix operator* (double s){
            return mult(s);
        }

        matrix flip(){
            matrix m;
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    m[i][j] = elems[n-i-1][n-j-1];
                }
            }
            return m;
        }
        
        int getTriangles(matrix* upperOut){
            matrix upper = copy();
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

                if (maxV < FLOAT_TOLERANCE) return -1; //Degenerate Matrix, Determinant 0

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
                    // lower[j][i] = factor;
                    for(int k = 0; k < n; k++){
                        upper[j][k] -= factor*upper[i][k];
                        
                    }
                }
            }
            *upperOut = upper;

            return numPivots;
        }

        double determinant(){
            matrix m;
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

        matrix getMinors(){
            double minorVals[n][n];
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
                    matrix<n-1,n-1> minor(vals);
                    minorVals[i][j] = minor.determinant();
                }
            }
            return matrix(minorVals);
        }

        matrix getCofactor(){
            matrix minors = getMinors();
            for(int i = 0; i < n; i++){
                for(int j = 0; j < n; j++){
                    if((i+j)%2){
                        minors[i][j]*=-1;
                    }
                }
            }
            return minors;
        }

        matrix getAdjugate(){
            return getCofactor().transpose();
        }

        matrix getInverse(){
            return getAdjugate() * (1/determinant());
        }
    };
}

namespace matlib{
    class vec3{
        public:
        double x,y,z;
        vec3(double x, double y, double z){
            this->x = x;
            this->y = y;
            this->z = z;
        }

        vec3(matrix<3,1> m){
            this->x = m[0][0];
            this->y = m[1][0];
            this->z = m[2][0];
        }

        matrix<3,1> getMatrix(){
            double p[3][1] = {{x},{y},{z}};
            return matrix<3,1>(p);
        }

        std::string to_string(){
            return '(' + std::to_string(x) + ',' + std::to_string(y) + ',' + std::to_string(z) + ')';
        }

        double dot(vec3 v){
            return x*v.x + y*v.y + z*v.z;
        }

        double mag(){
            return sqrt(dot(*this));
        }

        double angleBetween(vec3 v){
            return acos(dot(v)/mag()/v.mag());
        }

        vec3 cross(vec3 v){
            double x2 = y*v.z - z*v.y;
            double y2 = z*v.x - x*v.z;
            double z2 = x*v.y - y*v.x;
            return vec3(x2,y2,z2);
        }

        vec3 transform(matrix<3,3> &m){
            matrix<3,1> m2 = getMatrix();
            return vec3(m * m2);
        }
    };
    class vec4{
        public:
        double x,y,z,w;
        vec4(double x, double y, double z){
            this->x = x;
            this->y = y;
            this->z = z;
            this->w = 1;
        }

        vec4(vec3 v){
            x = v.x;
            y = v.y;
            z = v.z;
            w = 1;
        }

        vec4(matrix<4,1> m){
            this->x = m[0][0];
            this->y = m[1][0];
            this->z = m[2][0];
            this->w = m[3][0];
        }

        matrix<4,1> getMatrix(){
            double p[4][1] = {{x},{y},{z},{1}};
            return matrix<4,1>(p);
        }

        vec4 transform(matrix<4,4> &m){
            matrix<4,1> m2 = getMatrix();
            return vec4(m * m2);
        }
    };
}

namespace matlib{
    namespace conversion{
        double rad_to_deg(double r) {
            return r * (180.0 / M_PI);
        }

        double deg_to_rad(double d){
            return d * M_PI / 180.0;
        }
    }
}

namespace matlib{
    namespace transform{
        matrix<2,2> make2DRot(double a){
            double v[2][2] = {{cos(a),-sin(a)}, {sin(a), cos(a)}};
            return matrix<2,2>(v);
        }

        matrix<3,3> makeXRot(double a){
            double vals[3][3] = {{1,0,0},{0,cos(a),-sin(a)},{0,sin(a),cos(a)}};
            return matrix<3,3>(vals);
        }

        matrix<3,3> makeYRot(double a){
            double vals[3][3] = {{cos(a),0,sin(a)},{0,1,0},{-sin(a),0,cos(a)}};
            return matrix<3,3>(vals);
        }

        matrix<3,3> makeZRot(double a){
            double vals[3][3] = {{cos(a),-sin(a),0},{sin(a),cos(a),0},{0,0,1}};
            return matrix<3,3>(vals);
        }

        matrix<3,3> makeRot(double a, double b, double c){
            matrix<3,3> mx = makeXRot(a), my = makeYRot(b), mz = makeZRot(c);
            return mx * my * mz;
        }

        matrix<4,4> makeTranslation(vec3 v){
            double vals[4][4] = {{1,0,0,v.x},{0,1,0,v.y},{0,0,1,v.z},{0,0,0,1}};
            return matrix<4,4>(vals);
        }

        matrix<3,3> makeScale(double scale){
            return matrix<3,3>::Identity() * scale;
        }
    }
}

namespace matlib{
    template <>
    double matrix<0,0>::determinant(){
        return 0;
    }

    template <>
    double matrix<1,1>::determinant(){
        return elems[0][0];
    }

    template <>
    matrix<0,0> matrix<0,0>::getMinors(){
        return 0;
    }

    template <>
    matrix<1,1> matrix<1,1>::getMinors(){
        return 0;
    }

    template <>
    matrix<0,0> matrix<0,0>::getCofactor(){
        return 0;
    }

    template <>
    matrix<1,1> matrix<1,1>::getCofactor(){
        return 0;
    }

    template <>
    matrix<0,0> matrix<0,0>::getAdjugate(){
        return 0;
    }

    template <>
    matrix<1,1> matrix<1,1>::getAdjugate(){
        return 0;
    }

    template <>
    matrix<0,0> matrix<0,0>::getInverse(){
        return 0;
    }

    template <>
    matrix<1,1> matrix<1,1>::getInverse(){
        return 0;
    }

}