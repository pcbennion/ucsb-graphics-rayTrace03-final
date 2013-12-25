#ifndef _MATRIX3_H_
#define _MATRIX3_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "VEC3F.h"

using namespace std;

class MATRIX3 {
public:
  MATRIX3() {
    setToZero();
  };

  void setToZero() {
    for (int x = 0; x < 9; x++) 
      data[x] = 0;
  };

  void setToIdentity() {
    setToZero();
    data[0] = 1.0;
    data[4] = 1.0;
    data[8] = 1.0;
  };

  VEC3F operator*(const VEC3F& v) {
    VEC3F final;
    final.x = data[0] * v.x + data[1] * v.y + data[2] * v.z;
    final.y = data[3] * v.x + data[4] * v.y + data[5] * v.z;
    final.z = data[6] * v.x + data[7] * v.y + data[8] * v.z;
    return final;
  };

  static MATRIX3 I() {
    MATRIX3 final;
    final.setToIdentity();
    return final;
  };

  // note! C++ uses radians, but OpenGL uses degrees!
  static MATRIX3 rotateY(const float theta) {
    MATRIX3 final;
    final.data[0] = cos(theta);  final.data[1] = 0; final.data[2] = sin(theta);
    final.data[3] = 0;           final.data[4] = 1; final.data[5] = 0;
    final.data[6] = -sin(theta); final.data[7] = 0; final.data[8] = cos(theta);
    return final;
  };

  float data[9];
};
#endif
