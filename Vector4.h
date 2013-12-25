#ifndef VECTOR4_H
#define VECTOR4_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cassert>

class Vector4
{
 public:
  // Vector elements
  float elems[4];

  // Constructors
  Vector4()
    {elems[0]=0; elems[1]=0; elems[2]=0; elems[3]=0;};
  Vector4(float coords[4])
    {elems[0]=coords[0]; elems[1]=coords[1]; elems[2]=coords[2]; elems[3]=coords[3];};
  Vector4(float a, float b, float c, float d)
    {elems[0]=a; elems[1]=b; elems[2]=c; elems[3]=d;};
  Vector4(const Vector4& v)
    {elems[0]=v.elems[0]; elems[1]=v.elems[1]; elems[2]=v.elems[2]; elems[3]=v.elems[3];};

  // Vector operations:
  //		vector length, dot product, cross product
  void print();
  float dot(Vector4 v)
  {return (elems[0]*v[0] + elems[1]*v[1] + elems[2]*v[2] + elems[3]*v[3]);};

  // Subscript operator
  float& operator[] (const int i)
  {assert(i<4&&i>=0); return elems[i];};

  // Assignment operators
  Vector4 operator=(Vector4 rhs)
    {
      this->elems[0] = rhs[0]; this->elems[1] = rhs[1]; this->elems[2] = rhs[2]; this->elems[3] = rhs[3];
      return *this;
    }
  Vector4 operator=(float* rhs)
    {
      this->elems[0] = rhs[0]; this->elems[1] = rhs[1]; this->elems[2] = rhs[2]; this->elems[3] = rhs[3];
      return *this;
    }

  // Arithmetic operators:
  //		vector addition, subtraction, negation, and multiplying across elements
  //		scalar multiplication, division.
  Vector4 operator+=(Vector4 rhs)
  {elems[0]+=rhs[0]; elems[1]+=rhs[1]; elems[2]+=rhs[2]; elems[3]+=rhs[3]; return *this;};
  Vector4 operator+(Vector4 rhs)
  {return Vector4(elems[0] + rhs[0], elems[1] + rhs[1], elems[2] + rhs[2], elems[3] + rhs[3]);};
  Vector4 operator-()
  {elems[0]=-elems[0]; elems[1]=-elems[1]; elems[2]=-elems[2]; elems[3]=-elems[3]; return *this;};
  Vector4 operator-=(Vector4 rhs)
  {elems[0]-=rhs[0]; elems[1]-=rhs[1]; elems[2]-=rhs[2]; elems[3]-=rhs[3]; return *this;};
  Vector4 operator-(Vector4 rhs)
  {return Vector4(elems[0] - rhs[0], elems[1] - rhs[1], elems[2] - rhs[2], elems[3] - rhs[3]);};
  Vector4 operator*=(float rhs)
  {elems[0]*=rhs; elems[1]*=rhs; elems[2]*=rhs; elems[3]*=rhs; return *this;};
  Vector4 operator*(float rhs)
    {return Vector4(elems[0] * rhs, elems[1] * rhs, elems[2] * rhs, elems[3] * rhs);};
  Vector4 operator*(Vector4 rhs)
    {return Vector4(elems[0] * rhs[0], elems[1] * rhs[1], elems[2] * rhs[2], elems[3] * rhs[3]);};
  Vector4 operator/=(float rhs)
  {elems[0]/=rhs; elems[1]/=rhs; elems[2]/=rhs; elems[3]/=rhs; return *this;};
  Vector4 operator/(float rhs)
  {return Vector4(elems[0] / rhs, elems[1] / rhs, elems[2] / rhs, elems[3] / rhs);};

  // Comparison operators
  bool operator==(Vector4 rhs)
  {
    for(int i=0; i<4; i++){if(elems[i]!=rhs[i]) return false;}
    return true;
  };
  bool operator!=(Vector4 rhs)
  {
    for(int i=0; i<4; i++){if(elems[i]!=rhs[i]) return true;}
    return false;
  };
};

// Output operator for easy printing
static std::ostream& operator<<(std::ostream& os, Vector4 v)
{
  std::cout<<"<";
  for(int i=0; i<4; i++) {std::cout<<v[i]; if(i<3) std::cout<<", ";}
  std::cout<<">";
  return os;
}

#endif /* VECTOR4_H_ */
