#ifndef MATRIX4_H_
#define MATRIX4_H_

#include <math.h>
#include <iostream>
#include "VEC3F.h"
#include "VEC4F.h"

using namespace::std;

class MATRIX4
{
private:
  VEC4F row[4];
public:
  //
  // Initalization
  //
  MATRIX4() 
    { }
  MATRIX4(float* data)
    {
      for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
          row[i][j] = data[i * 4 + j];
    }
  MATRIX4(double* data)
    {
      for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
          row[i][j] = data[i * 4 + j];
    }
  MATRIX4(const MATRIX4 &M) 
    { *this = M; }
  MATRIX4(const VEC3F &r0, const VEC3F &r1, const VEC3F &r2)
    { 
      row[0] = VEC4F(r0, true);
      row[1] = VEC4F(r1, true);
      row[2] = VEC4F(r2, true);
      row[3][3] = 1.0; 
    }
  MATRIX4(const VEC4F &r0, const VEC4F &r1, const VEC4F &r2, const VEC4F &r3)
    { row[0] = r0, row[1] = r1, row[2] = r2, row[3] = r3; }
  //
  // Assignment
  //
  MATRIX4& operator = (const MATRIX4 &M) 
    { row[0] = M[0]; row[1] = M[1]; row[2] = M[2]; row[3] = M[3]; return *this; }


  MATRIX4 transpose() const 
    { 
      return MATRIX4(col(0), col(1), col(2), col(3)); 
    } 

  //
  // Vector transforms
  //
  inline VEC3F transformPoint(const VEC3F &point) const;

  inline VEC3F transformDirection(const VEC3F &direction) const;

  //
  // Accessors
  //
  inline VEC4F& operator[] (int r)                 { return row[r]; }
  inline const VEC4F& operator[] (int r) const     { return row[r]; }
  inline operator float*()                          { return row[0]; }
  inline operator const float*()                    { return row[0]; }
  inline operator const float*() const              { return row[0]; }
  float& operator()(const int x, const int y)  { return row[y][x]; }
  
  inline void setColumn(int c, const VEC4F &values)
    {
      row[0][c] = values[0];
      row[1][c] = values[1];
      row[2][c] = values[2];
      row[3][c] = values[3];
    }
  inline void setRow(int r, const VEC4F &values)
    { row[r] = values; }
  inline VEC4F col(int c) const
    {
      return VEC4F(row[0][c], row[1][c], row[2][c], row[3][c]);
    }

  //
  // Transformation matrices
  //
  static MATRIX4 I();
  static MATRIX4 Scaling(const VEC3F &s);
  static MATRIX4 Scaling(float s);
  static MATRIX4 Translation(const VEC3F &pos);
  static MATRIX4 Rotation(const VEC3F &axis, float angle);
  static MATRIX4 Rotation(const VEC3F &basis1, const VEC3F &basis2, const VEC3F &basis3);
  static MATRIX4 RotationX(float theta);
  static MATRIX4 RotationY(float theta);
  static MATRIX4 RotationZ(float theta);
  static MATRIX4 LookAt(const VEC3F &Eye, const VEC3F &At, const VEC3F &Up);
  static MATRIX4 Perspective(float FOV, float Aspect, float ZNear, float ZFar);
};

//
// Overloaded operators
//
inline VEC4F operator*(const MATRIX4 &Left, const VEC4F &Right)
{
  return VEC4F(Left[0] * Right, Left[1] * Right, Left[2] * Right, Left[3] * Right);
}

inline ostream& operator << (ostream &os, const MATRIX4 &m)
{
  return os << "[" << endl << m[0] << endl << m[1] << endl << m[2] << endl << m[3] << "]" << endl;
}
inline VEC3F MATRIX4::transformPoint(const VEC3F &point) const
{
  VEC4F hPoint(point, false);
  return ((*this) * hPoint).toVEC3F(true);
}

inline VEC3F MATRIX4::transformDirection(const VEC3F &direction) const
{
  VEC4F hDirection(direction, true);
  return ((*this) * hDirection).toVEC3F(false);
}

MATRIX4 operator * (const MATRIX4 &Left, const MATRIX4 &Right);
MATRIX4 operator + (const MATRIX4 &Left, const MATRIX4 &Right);
MATRIX4 operator - (const MATRIX4 &Left, const MATRIX4 &Right);
#endif
