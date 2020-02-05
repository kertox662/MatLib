#ifndef MATLIBTRANSFORM_H //MATLIBTRANSFORM_H Header Guard
#define MATLIBTRANSFORM_H

#include "matlibMat.hpp"
#include "matlibVec.hpp"

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

#endif //MATLIBTRANSFORM_H Header Guard