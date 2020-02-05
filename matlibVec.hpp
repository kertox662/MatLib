#ifndef MATLIBVEC_H //MATLIBVEC_H Header Guard
#define MATLIBVEC_H

#include "matlibMat.hpp"

namespace Matlib{
    class vec3;
    class vec4;
}

namespace Matlib{ //Vector Class Member Declaration
    class vec3{
        public:
        double x,y,z;
        
        vec3();
        vec3(double x, double y, double z);
        vec3(double v[3]);
        vec3(Matrix<3,1> m);
        vec3(const vec3 &v);

        Matrix<3,1> getMatrix();
        std::string to_string();

        vec3 add(vec3 v);
        vec3 sub(vec3 v);
        vec3 neg();
        vec3 operator+(vec3 v);
        vec3 operator-(vec3);
        vec3 operator-();

        vec3 mult(double s);
        vec3 operator*(double s);

        double mag();

        double dot(vec3 v);
        vec3 cross(vec3 v);

        double angleBetween(vec3 v);

        vec3 projectOnto(vec3);

        vec3 transform(Matrix<3,3> &m);
        vec3 transform(Matrix<4,4> &m);
    };
    class vec4{
        public:
        double x,y,z,w;
        vec4();
        vec4(double x, double y, double z);
        vec4(double x, double y, double z, double w);
        vec4(double v[4]);
        vec4(vec3 v);
        vec4(const vec4 &v);
        vec4(Matrix<4,1> m);

        vec4 add(vec4 v);
        vec4 sub(vec4 v);
        vec4 mult(double s);
        vec4 neg();
        vec4 operator+(vec4 v);
        vec4 operator-(vec4);
        vec4 operator-();
        vec4 operator*(double s);

        double mag();

        double dot(vec4 v);

        Matrix<4,1> getMatrix();
        std::string to_string();

        vec4 transform(Matrix<4,4> &m);
    };
}

namespace Matlib{ //Vector 3D Class Member Definitions
    vec3::vec3(){
        x = y = z = 0;
    }

    vec3::vec3(double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    vec3::vec3(double v[3]){
        x = v[0];
        y = v[1];
        z = v[2];
    }

    vec3::vec3(Matrix<3,1> m){
        this->x = m[0][0];
        this->y = m[1][0];
        this->z = m[2][0];
    }

    vec3::vec3(const vec3 &v){
        x = v.x;
        y = v.y;
        z = v.z;
    }

    Matrix<3,1> vec3::getMatrix(){
        double p[3][1] = {{x},{y},{z}};
        return Matrix<3,1>(p);
    }

    std::string vec3::to_string(){
        return '(' + std::to_string(x) + ',' + std::to_string(y) + ',' + std::to_string(z) + ')';
    }

    vec3 vec3::add(vec3 v){
        vec3 sum;
        sum.x = x + v.x;
        sum.y = y + v.y;
        sum.z = z + v.y;
        return sum;
    }

    vec3 vec3::operator+(vec3 v){
        return add(v);
    }

    vec3 vec3::neg(){
        vec3 v;
        v.x = -x;
        v.y = -y;
        v.z = -z;
        return v;
    }

    vec3 vec3::operator-(){
        return this->neg();
    }

    vec3 vec3::sub(vec3 v){
        return this->add(-v);
    }

    vec3 vec3::operator-(vec3 v){
        return sub(v);
    }

    vec3 vec3::mult(double s){
        vec3 v(*this);
        v.x *= s;
        v.y *= s;
        v.z *= s;
        return v;
    }

    vec3 vec3::operator*(double s){
        return mult(s);
    }

    double vec3::dot(vec3 v){
        return x*v.x + y*v.y + z*v.z;
    }

    double vec3::mag(){
        return sqrt(dot(*this));
    }

    double vec3::angleBetween(vec3 v){
        return acos(dot(v)/mag()/v.mag());
    }

    vec3 vec3::cross(vec3 v){
        double x2 = y*v.z - z*v.y;
        double y2 = z*v.x - x*v.z;
        double z2 = x*v.y - y*v.x;
        return vec3(x2,y2,z2);
    }

    vec3 vec3::transform(Matrix<3,3> &m){
        Matrix<3,1> m2 = getMatrix();
        return vec3(m * m2);
    }

    vec3 vec3::transform(Matrix<4,4> &m){
        vec4 v(*this);
        vec4 v2 = v.transform(m);
        vec3 result;
        result.x = v.x;
        result.y = v.y;
        result.z = v.y;
        return result;
    }

    vec3 vec3::projectOnto(vec3 v){
        double s = dot(v) / v.dot(v);
        return v*s;
    }
}
namespace Matlib{ //Vector 4D Class Member Definitions
    vec4::vec4(){
        x = y = z = w = 0;
    }
    
    vec4::vec4(const vec4 &v){
        x = v.x;
        y = v.y;
        z = v.z;
        w = v.w;
    }

    vec4::vec4(double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = 1;
    }

    vec4::vec4(double x, double y, double z, double w){
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = w;
    }

    vec4::vec4(double v[4]){
        x = v[0];
        y = v[1];
        z = v[2];
        w = v[3];
    }

    vec4::vec4(vec3 v){
        x = v.x;
        y = v.y;
        z = v.z;
        w = 1;
    }

    vec4::vec4(Matrix<4,1> m){
        this->x = m[0][0];
        this->y = m[1][0];
        this->z = m[2][0];
        this->w = m[3][0];
    }

    Matrix<4,1> vec4::getMatrix(){
        double p[4][1] = {{x},{y},{z},{1}};
        return Matrix<4,1>(p);
    }

    std::string vec4::to_string(){
        return "(" + std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z) + "," + std::to_string(w) + ")";
    }

    vec4 vec4::add(vec4 v){
        vec4 result;
        result.x = x + v.x;
        result.y = x + v.y;
        result.z = x + v.z;
        result.w = x + v.w;
        return result;
    }

    vec4 vec4::operator+(vec4 v){
        return this->add(v);
    }

    vec4 vec4::neg(){
        vec4 v;
        v.x = -x;
        v.y = -y;
        v.z = -z;
        v.w = -w;
        return v;
    }

    vec4 vec4::operator-(){
        return this->neg();
    }

    vec4 vec4::sub(vec4 v){
        return this->add(-v);
    }
    
    vec4 vec4::operator-(vec4 v){
        return this->sub(v);
    }

    vec4 vec4::mult(double s){
        vec4 v(*this);
        v.x *= s;
        v.y *= s;
        v.z *= s;
        v.w *= s;
        return v;
    }

    vec4 vec4::operator*(double s){
        return this->mult(s);
    }

    double vec4::mag(){
        return sqrt(x*x + y*y + z*z + w*w);
    }

    double vec4::dot(vec4 v){
        double s = 0;
        s += x * v.x;
        s += y * v.y;
        s += z * v.z;
        s += w * v.w;
        return s;
    }

    vec4 vec4::transform(Matrix<4,4> &m){
        Matrix<4,1> m2 = getMatrix();
        return vec4(m * m2);
    }
}

#endif //MATLIBVEC_H Header Guard