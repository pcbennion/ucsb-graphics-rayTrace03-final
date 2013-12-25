#ifndef VEC3F_H
#define VEC3F_H

#include <iostream>
#include <math.h>

 
class VEC3F{
public:
  union {
    struct { float x,y,z; };
    float element[3];
  };
  VEC3F(float s = 0.f) { element[0] = element[1] = element[2] = s; }
  VEC3F(float x, float y, float z) {
    element[0] = x; element[1] = y; element[2] = z;
  }

  // Copy constructors & assignment operators
  VEC3F(const VEC3F& v) { *this = v; }
  VEC3F(const float v[3])
    { element[0] = v[0]; element[1] = v[1]; element[2] = v[2];}

  VEC3F& operator=(const VEC3F& v)
    { element[0] = v[0]; element[1] = v[1]; element[2] = v[2]; return *this; }
  VEC3F& operator=(float s) { element[0] = element[1] = element[2] = s; return *this; }

  // Access methods
  operator       float*()       { return element; }
  operator const float*() const { return element; }

  float& operator[](int i)       { return element[i]; }
  float  operator[](int i) const { return element[i]; }

  // Assignment and in-place arithmetic methods
  inline VEC3F& operator+=(const VEC3F& v);
  inline VEC3F& operator-=(const VEC3F& v);
  inline VEC3F& operator*=(float s);
  inline VEC3F& operator/=(float s);

  void normalize() {
    float l = element[0] * element[0] + element[1] * element[1] + element[2] * element[2];
    if( l!=1.0 && l!=0.0 )  *this /= sqrt(l);
  };
};
////////////////////////////////////////////////////////////////////////
// Method definitions
////////////////////////////////////////////////////////////////////////

inline VEC3F& VEC3F::operator+=(const VEC3F& v)
  { element[0] += v[0];   element[1] += v[1];   element[2] += v[2]; return *this; }

inline VEC3F& VEC3F::operator-=(const VEC3F& v)
  { element[0] -= v[0];   element[1] -= v[1];   element[2] -= v[2]; return *this; }

inline VEC3F& VEC3F::operator*=(float s)
  { element[0] *= s;   element[1] *= s;   element[2] *= s; return *this; }

inline VEC3F& VEC3F::operator/=(float s)
  { element[0] /= s;   element[1] /= s;   element[2] /= s; return *this; }


////////////////////////////////////////////////////////////////////////
// Operator definitions
////////////////////////////////////////////////////////////////////////

inline VEC3F operator+(const VEC3F &u, const VEC3F& v)
  { return VEC3F(u[0]+v[0], u[1]+v[1], u[2]+v[2]); }


inline VEC3F operator-(const VEC3F &u, const VEC3F& v)
  { return VEC3F(u[0]-v[0], u[1]-v[1], u[2]-v[2]); }

inline VEC3F operator-(const VEC3F &v)
  { return VEC3F(-v[0], -v[1], -v[2]); }

inline VEC3F operator*(float s, const VEC3F &v)
  { return VEC3F(v[0]*s, v[1]*s, v[2]*s); }
inline VEC3F operator*(const VEC3F &v, float s)
  { return s*v; }

inline VEC3F operator/(const VEC3F &v, float s)
  { return VEC3F(v[0]/s, v[1]/s, v[2]/s); }

inline float operator*(const VEC3F &u, const VEC3F& v)
  { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }

inline VEC3F cross(const VEC3F& u, const VEC3F& v){ 
  return VEC3F(u[1]*v[2] - v[1]*u[2], 
              -u[0]*v[2] + v[0]*u[2], 
               u[0]*v[1] - v[0]*u[1]); 
}

inline std::ostream &operator<<(std::ostream &out, const VEC3F& v)
  { return out << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")"; }

////////////////////////////////////////////////////////////////////////
// Misc. function definitions
////////////////////////////////////////////////////////////////////////

inline float norm2(const VEC3F& v)  { return v*v; }
inline float norm(const VEC3F& v)   { return sqrt(norm2(v)); }
inline VEC3F normalize(const VEC3F& v) 
{
  float l = norm(v);
  if(l != 0.f)
    return VEC3F(v[0] / l, v[1] / l, v[2] / l);
  else
    return VEC3F(0.f);
}
#endif
