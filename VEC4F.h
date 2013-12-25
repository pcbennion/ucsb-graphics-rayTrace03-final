#ifndef VEC4F_H
#define VEC4F_H

#include <iostream>
#include <math.h>
#include "VEC3F.h"

class VEC4F{
private:
  float element[4];
public:
  VEC4F(float s = 0) { element[0] = element[1] = element[2] = element[3] = s; }
  VEC4F(float x, float y, float z, float a = 1) {
    element[0] = x; element[1] = y; element[2] = z; element[3] = a;
  }

  // Copy constructors & assignment operators
  VEC4F(const VEC3F& v, bool isDirection){ 
    element[0] = v[0]; element[1] = v[1]; element[2] = v[2]; 
    if(isDirection)
      element[3] = 0; 
    else
      element[3] = 1;
  }
  VEC4F(const VEC4F& v) { *this = v; }
  VEC4F(const float v[4])
    { element[0] = v[0]; element[1] = v[1]; element[2] = v[2]; element[3] = v[3]; }
  VEC4F& operator=(const VEC4F& v)
    { element[0] = v[0]; element[1] = v[1]; element[2] = v[2]; element[3] = v[3]; return *this; }

  // Access methods
  operator       float*()       { return element; }
  operator const float*() const { return element; }

  float& operator[](int i)       { return element[i]; }
  float  operator[](int i) const { return element[i]; }

  inline VEC3F toVEC3F(bool regularize) const{ 
    if(!regularize)
      return VEC3F(element[0], element[1], element[2]); 
    else if(fabs(element[3]) != 0)
      return VEC3F(element[0] / element[3], element[1] / element[3], element[2] / element[3]); 
    else 
      return VEC3F(0.f);
  }
};

inline float operator*(const VEC4F &u, const VEC4F& v)
  { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2] + u[3]*v[3]; }
inline VEC4F operator+(const VEC4F &u, const VEC4F& v)
  { return VEC4F(u[0]+v[0], u[1]+v[1], u[2]+v[2], u[3]+v[3]); }
inline VEC4F operator-(const VEC4F &u, const VEC4F& v)
  { return VEC4F(u[0]-v[0], u[1]-v[1], u[2]-v[2], u[3]-v[3]); }

inline std::ostream &operator<<(std::ostream &out, const VEC4F& v)
  { return out << "(" << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << ")"; }

#endif
