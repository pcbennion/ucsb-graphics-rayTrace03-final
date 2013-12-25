#ifndef MATRIX3_H
#define MATRIX3_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cassert>

#include "Vector3.h"
#include "MATRIX3.h"
#include "MATRIX4.h"

class Matrix3
{
public:
    // Matrix elements
    float elems[9];

    // Constructors
    Matrix3()
        {for(int i=0;i<9;i++) elems[i]=0;};
    Matrix3(float a, float b, float c, float d, float e, float f, float g, float h, float i)
    {
        elems[0]=a;elems[1]=b;elems[2]=c;
        elems[3]=d;elems[4]=e;elems[5]=f;
        elems[6]=g;elems[7]=h;elems[8]=i;
    }
    Matrix3(Vector3 x, Vector3 y, Vector3 z)
        {for(int i=0;i<3;i++) {elems[i]=x[i]; elems[i+3]=y[i]; elems[i+6]=z[i];}};

    Matrix3(float* mat)
        {for(int i=0;i<9;i++) elems[i]=mat[i];};

    // Transpose function
    Matrix3 transpose()
    {
        Matrix3 out = Matrix3();
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
                out[3*j+i] = this->elems[3*i+j];
        }
        return out;
    };

    // Subscript operator
    float& operator[] (const int i)
        {assert(i<9&&i>=0); return elems[i];};

    // Assignment operators
    Matrix3 operator=(Matrix3 rhs)
    {
        for(int i=0;i<9;i++) elems[i]=rhs[i];
        return *this;
    }
    Matrix3 operator=(float* rhs)
    {
        for(int i=0;i<9;i++) elems[i]=rhs[i];
        return *this;
    }
    Matrix3 operator=(MATRIX3 rhs)
    {
        for(int i=0;i<9;i++) elems[i]=rhs.data[i];
        return *this;
    }
    Matrix3 operator=(MATRIX4 rhs)
    {
        for(int i=0;i<9;i++) elems[i]=rhs[i/3][i%3];
        return *this;
    }

    // Simple Arithmetic operators:
    //		matrix addition, subtraction, negation
    //		scalar multiplication
    Matrix3 operator+=(Matrix3 rhs)
        {for(int i=0;i<9;i++) elems[i]+=rhs[i]; return *this;};
    Matrix3 operator+(Matrix3 rhs)
    {
        Matrix3 out = *this;
        out += rhs;
        return out;
    };
    Matrix3 operator-=(Matrix3 rhs)
        {for(int i=0;i<9;i++) elems[i]-=rhs[i]; return *this;};
    Matrix3 operator-(Matrix3 rhs)
    {
        Matrix3 out = *this;
        out -= rhs;
        return out;
    };
    Matrix3 operator*=(float rhs)
        {for(int i=0;i<9;i++) elems[i]*=rhs; return *this;};
    Matrix3 operator*(float rhs)
    {
        Matrix3 out = *this;
        out *= rhs;
        return out;
    };

    // More Arithmetic Operators:
    //      multiplication w/ matricies, vectors
    Matrix3 operator*(Matrix3 rhs)
    {
        Matrix3 out = Matrix3();
        int index = 0;
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                out[index] = this->elems[3*i+0]*rhs[3*0+j]
                            +this->elems[3*i+1]*rhs[3*1+j]
                            +this->elems[3*i+2]*rhs[3*2+j];
                index++;
            }
        }
        return out;
    };
    Vector3 operator*(Vector3 rhs)
    {
        Vector3 out = Vector3();
        int index = 0;
        for(int i=0; i<3; i++)
        {
            out[index] = this->elems[3*i+0]*rhs[0]
                        +this->elems[3*i+1]*rhs[1]
                        +this->elems[3*i+2]*rhs[2];
            index++;
        }
        return out;
    };

    // Comparison operators
    bool operator==(Matrix3 rhs)
    {
        for(int i=0; i<9; i++){if(elems[i]!=rhs[i]) return false;}
        return true;
    };
    bool operator!=(Matrix3 rhs)
    {
        for(int i=0; i<9; i++){if(elems[i]!=rhs[i]) return true;}
        return false;
    };
};

// Output operator for easy printing
static std::ostream& operator<<(std::ostream& os, Matrix3 v)
{
    std::cout<<std::endl;
    for(int i=0; i<9; i++)
    {
        if(i%3==0) std::cout<<"[";
        std::cout<<v[i];
        if((i+1)%3==0) {std::cout<<"]";std::cout<<std::endl;}
        else std::cout<<", ";
    }
    std::cout<<std::endl;
    return os;
}

#endif // MATRIX3_H
