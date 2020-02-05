#ifndef MATLIBMATM_H //MATLIBMATM_H Header Guard
#define MATLIBMATM_H

#include <string>
#include <cmath>

#ifndef MATLIB_FLOAT_TOLERANCE
#define MATLIB_FLOAT_TOLERANCE 0.000000001
#endif

namespace Matlib{ //Class Declaration
    class MatrixM;
}

namespace Matlib{ //Class Member Declaration
    class MatrixM{
        public:
        struct Size{
            int rows,columns;
            std::string to_string(){
                return std::to_string(rows) + 'x' + std::to_string(columns);
            }
        };

        private:
        Size dimensions;
        double* elems;

        public:
        MatrixM(int r, int c);
        MatrixM(int r, int c, double* e);
        MatrixM(const MatrixM &m);
        MatrixM operator=(const MatrixM &m);

        MatrixM IdentityLeft();
        MatrixM IdentityRight();

        double* operator[](int n);

        std::string to_string();

        MatrixM copy();

        MatrixM transpose();

        MatrixM add(MatrixM &m);
        MatrixM sub(MatrixM &m);
        MatrixM neg();
        MatrixM mult(MatrixM &m);
        MatrixM mult(double s);

        MatrixM operator+(MatrixM &m);
        MatrixM operator-(MatrixM &m);
        MatrixM operator-();
        MatrixM operator*(MatrixM &m);
        MatrixM operator*(double s);

        double determinant();

        MatrixM getMinors();
        MatrixM getCofactors();
        MatrixM getAdjugate();
        MatrixM Inverse();
        MatrixM Inverse(int* err);

        private:
        int getTriangles(MatrixM* upperOut){
            if(dimensions.rows != dimensions.columns) return -1;
            int n = dimensions.rows;
            
            MatrixM upper = copy();
            int P[n]; //Permutation Matrix;
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
                    int permSwap = P[i]; //Permuation Matrix swapper
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

namespace Matlib{ //Class Member Definitions
    MatrixM::MatrixM(int r, int c){
        dimensions = {r,c};
    }

    MatrixM::MatrixM(int r, int c, double* e){
        dimensions = {r,c};
        elems = e;
    }
}

#endif //MATLIBMATM_H Header Guard