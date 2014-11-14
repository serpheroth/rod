#ifndef SVD_H
#define SVD_H
#include <iostream>
#include <math.h>

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

#include <math.h>
#include "MY_MATH.h"
#include "vector_lib.h"
#define NR_SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
namespace git {

#define _gamma 5.828427124f // FOUR_GAMMA_SQUARED = sqrt(8)+3;
#define _cstar 0.923879532f // cos(pi/8)
#define _sstar 0.3826834323f // sin(p/8)
const float kEpsilon = 1e-6f;

/* This is a novel and fast routine for the reciprocal square root of an
IEEE float (single precision).

http://www.lomont.org/Math/Papers/2003/InvSqrt.pdf
http://playstation2-linux.com/download/p2lsd/fastrsqrt.pdf
http://www.beyond3d.com/content/articles/8/
*/
inline float rsqrt(float x) {
  return 1.0f / sqrtf(x);
  // int ihalf = *(int *)&x - 0x00800000; // Alternative to next line,
  // float xhalf = *(float *)&ihalf;      // for sufficiently large nos.
  float xhalf = 0.5f * x;
  //  int i = *(int *)&x;          // View x as an int.
  static_assert(sizeof(int) == sizeof(float), "int and float do not have same size");
  int i = 0;
  i = *((int*) &x);
  // i = 0x5f3759df - (i >> 1);   // Initial guess (traditional).
  i = 0x5f375a82 - (i >> 1);   // Initial guess (slightly better).
  //  x = *(float *)&i;            // View i as float.
  x = *((float*) &i);
//  memcpy(&x, &i, sizeof(float));
  x = x * (1.5f - xhalf * x * x); // Newton step.
  // x = x*(1.5008908 - xhalf*x*x);  // Newton step for a balanced error.
  return x;
}

/* This is rsqrt with an additional step of the Newton iteration, for
increased accuracy. The constant 0x5f37599e makes the relative error
range from 0 to -0.00000463.
   You can't balance the error by adjusting the constant. */
inline float rsqrt1(float x) {
  return 1.0f / sqrtf(x);
  float xhalf = 0.5f * x;
  //  int i = *(int *)&x;          // View x as an int.
  int i = 0;
//  memcpy(&i, &x, sizeof(int));
  i = *((int*) &x);
  i = 0x5f37599e - (i >> 1);   // Initial guess.
  //  x = *(float *)&i;            // View i as float.
//  memcpy(&x, &i, sizeof(float));
  x = *((float*) &i);
  x = x * (1.5f - xhalf * x * x); // Newton step.
  x = x * (1.5f - xhalf * x * x); // Newton step again.
  return x;
}

inline float accurateSqrt(float x) {
  return sqrtf(x);
  return x * rsqrt1(x);
}

inline void condSwap(bool c, float &X, float &Y) {
  // used in step 2
  float Z = X;
  X = c ? Y : X;
  Y = c ? Z : Y;
}

inline void condNegSwap(bool c, float &X, float &Y) {
  // used in step 2 and 3
  float Z = -X;
  X = c ? Y : X;
  Y = c ? Z : Y;
}

// matrix multiplication M = A * B
inline void multAB(const float* a, const float* b, float* m) {
  m[0] = a[0] * b[0] + a[1] * b[3] + a[2] * b[6];
  m[1] = a[0] * b[1] + a[1] * b[4] + a[2] * b[7];
  m[2] = a[0] * b[2] + a[1] * b[5] + a[2] * b[8];

  m[3] = a[3] * b[0] + a[4] * b[3] + a[5] * b[6];
  m[4] = a[3] * b[1] + a[4] * b[4] + a[5] * b[7];
  m[5] = a[3] * b[2] + a[4] * b[5] + a[5] * b[8];

  m[6] = a[6] * b[0] + a[7] * b[3] + a[8] * b[6];
  m[7] = a[6] * b[1] + a[7] * b[4] + a[8] * b[7];
  m[8] = a[6] * b[2] + a[7] * b[5] + a[8] * b[8];
}

// matrix multiplication M = Transpose[A] * B
inline void multAtB(const float* a, const float* b, float* m) {
  m[0] = a[0] * b[0] + a[3] * b[3] + a[6] * b[6];
  m[1] = a[0] * b[1] + a[3] * b[4] + a[6] * b[7];
  m[2] = a[0] * b[2] + a[3] * b[5] + a[6] * b[8];

  m[3] = a[1] * b[0] + a[4] * b[3] + a[7] * b[6];
  m[4] = a[1] * b[1] + a[4] * b[4] + a[7] * b[7];
  m[5] = a[1] * b[2] + a[4] * b[5] + a[7] * b[8];

  m[6] = a[2] * b[0] + a[5] * b[3] + a[8] * b[6];
  m[7] = a[2] * b[1] + a[5] * b[4] + a[8] * b[7];
  m[8] = a[2] * b[2] + a[5] * b[5] + a[8] * b[8];
}

inline void quatToMat3(const float * qV,
                       float* m
                      ) {
  float w = qV[3];
  float x = qV[0];
  float y = qV[1];
  float z = qV[2];

  float qxx = x * x;
  float qyy = y * y;
  float qzz = z * z;
  float qxz = x * z;
  float qxy = x * y;
  float qyz = y * z;
  float qwx = w * x;
  float qwy = w * y;
  float qwz = w * z;

  m[0] = 1 - 2 * (qyy + qzz);
  m[1] = 2 * (qxy - qwz);
  m[2] = 2 * (qxz + qwy);
  m[3] = 2 * (qxy + qwz);
  m[4] = 1 - 2 * (qxx + qzz);
  m[5] = 2 * (qyz - qwx);
  m[6] = 2 * (qxz - qwy);
  m[7] = 2 * (qyz + qwx);
  m[8] = 1 - 2 * (qxx + qyy);
}

inline void approximateGivensQuaternion(float a11, float a12, float a22, float &ch, float &sh) {
  /*
       * Given givens angle computed by approximateGivensAngles,
       * compute the corresponding rotation quaternion.
       */
  ch = 2 * (a11 - a22);
  sh = a12;
  bool b = _gamma * sh * sh < ch * ch;
  // fast rsqrt function suffices
  // rsqrt2 (https://code.google.com/p/lppython/source/browse/algorithm/HDcode/newCode/rsqrt.c?r=26)
  // is even faster but results in too much error
  float w = rsqrt(ch * ch + sh * sh);
  ch = b ? w * ch : _cstar;
  sh = b ? w * sh : _sstar;
}

inline void jacobiConjugation( const int x, const int y, const int z,
                               float &s11,
                               float &s21, float &s22,
                               float &s31, float &s32, float &s33,
                               float * qV) {
  float ch, sh;
  approximateGivensQuaternion(s11, s21, s22, ch, sh);

  float scale = ch * ch + sh * sh;
  float a = (ch * ch - sh * sh) / scale;
  float b = (2 * sh * ch) / scale;

  // make temp copy of S
  float _s11 = s11;
  float _s21 = s21; float _s22 = s22;
  float _s31 = s31; float _s32 = s32; float _s33 = s33;

  // perform conjugation S = Q'*S*Q
  // Q already implicitly solved from a, b
  s11 = a * (a * _s11 + b * _s21) + b * (a * _s21 + b * _s22);
  s21 = a * (-b * _s11 + a * _s21) + b * (-b * _s21 + a * _s22);	s22 = -b * (-b * _s11 + a * _s21) + a * (-b * _s21 + a * _s22);
  s31 = a * _s31 + b * _s32;								s32 = -b * _s31 + a * _s32; s33 = _s33;

  // update cumulative rotation qV
  float tmp[3];
  tmp[0] = qV[0] * sh;
  tmp[1] = qV[1] * sh;
  tmp[2] = qV[2] * sh;
  sh *= qV[3];

  qV[0] *= ch;
  qV[1] *= ch;
  qV[2] *= ch;
  qV[3] *= ch;

  // (x,y,z) corresponds to ((0,1,2),(1,2,0),(2,0,1))
  // for (p,q) = ((0,1),(1,2),(0,2))
  qV[z] += sh;
  qV[3] -= tmp[z]; // w
  qV[x] += tmp[y];
  qV[y] -= tmp[x];

  // re-arrange matrix for next iteration
  _s11 = s22;
  _s21 = s32; _s22 = s33;
  _s31 = s21; _s32 = s31; _s33 = s11;
  s11 = _s11;
  s21 = _s21; s22 = _s22;
  s31 = _s31; s32 = _s32; s33 = _s33;

}

inline float dist2(float x, float y, float z) {
  return x * x + y * y + z * z;
}

// finds transformation that diagonalizes a symmetric matrix
inline void jacobiEigenanlysis( // symmetric matrix
  float &s11,
  float &s21, float &s22,
  float &s31, float &s32, float &s33,
  // quaternion representation of V
  float * qV) {
  qV[3] = 1; qV[0] = 0; qV[1] = 0; qV[2] = 0; // follow same indexing convention as GLM
  for (int i = 0; i < 4; i++) {
    // we wish to eliminate the maximum off-diagonal element
    // on every iteration, but cycling over all 3 possible rotations
    // in fixed order (p,q) = (1,2) , (2,3), (1,3) still retains
    //  asymptotic convergence
    jacobiConjugation(0, 1, 2, s11, s21, s22, s31, s32, s33, qV); // p,q = 0,1
    jacobiConjugation(1, 2, 0, s11, s21, s22, s31, s32, s33, qV); // p,q = 1,2
    jacobiConjugation(2, 0, 1, s11, s21, s22, s31, s32, s33, qV); // p,q = 0,2
  }
}


inline void sortSingularValues(// matrix that we want to decompose
  float &b11, float &b12, float &b13,
  float &b21, float &b22, float &b23,
  float &b31, float &b32, float &b33,
  // sort V simultaneously
  float &v11, float &v12, float &v13,
  float &v21, float &v22, float &v23,
  float &v31, float &v32, float &v33) {
  float rho1 = dist2(b11, b21, b31);
  float rho2 = dist2(b12, b22, b23);
  float rho3 = dist2(b13, b23, b33);
  bool c;
  c = rho1 < rho2;
  condNegSwap(c, b11, b12); condNegSwap(c, v11, v12);
  condNegSwap(c, b21, b22); condNegSwap(c, v21, v22);
  condNegSwap(c, b31, b32); condNegSwap(c, v31, v32);
  condSwap(c, rho1, rho2);
  c = rho1 < rho3;
  condNegSwap(c, b11, b13); condNegSwap(c, v11, v13);
  condNegSwap(c, b21, b23); condNegSwap(c, v21, v23);
  condNegSwap(c, b31, b33); condNegSwap(c, v31, v33);
  condSwap(c, rho1, rho3);
  c = rho2 < rho3;
  condNegSwap(c, b12, b13); condNegSwap(c, v12, v13);
  condNegSwap(c, b22, b23); condNegSwap(c, v22, v23);
  condNegSwap(c, b32, b33); condNegSwap(c, v32, v33);
}

inline
void QRGivensQuaternion(float a1, float a2, float &ch, float &sh) {
  // a1 = pivot point on diagonal
  // a2 = lower triangular entry we want to annihilate
  float epsilon = kEpsilon;//EPSILON;
  float rho = accurateSqrt(a1 * a1 + a2 * a2);

  sh = rho > epsilon ? a2 : 0;
  ch = fabs(a1) + dj::Max(rho, epsilon);
  bool b = a1 < 0;
  condSwap(b, sh, ch);
  float w = rsqrt(ch * ch + sh * sh);
  ch *= w;
  sh *= w;
}


inline void QRDecomposition(// matrix that we want to decompose
  float b11, float b12, float b13,
  float b21, float b22, float b23,
  float b31, float b32, float b33,
  // output Q
  float &q11, float &q12, float &q13,
  float &q21, float &q22, float &q23,
  float &q31, float &q32, float &q33,
  // output R
  float &r11, float &r12, float &r13,
  float &r21, float &r22, float &r23,
  float &r31, float &r32, float &r33) {
  float ch1, sh1, ch2, sh2, ch3, sh3;
  float a, b;

  // first givens rotation (ch,0,0,sh)
  QRGivensQuaternion(b11, b21, ch1, sh1);
  a = 1 - 2 * sh1 * sh1;
  b = 2 * ch1 * sh1;
  // apply B = Q' * B
  r11 = a * b11 + b * b21;  r12 = a * b12 + b * b22;  r13 = a * b13 + b * b23;
  r21 = -b * b11 + a * b21; r22 = -b * b12 + a * b22; r23 = -b * b13 + a * b23;
  r31 = b31;          r32 = b32;          r33 = b33;

  // second givens rotation (ch,0,-sh,0)
  QRGivensQuaternion(r11, r31, ch2, sh2);
  a = 1 - 2 * sh2 * sh2;
  b = 2 * ch2 * sh2;
  // apply B = Q' * B;
  b11 = a * r11 + b * r31;  b12 = a * r12 + b * r32;  b13 = a * r13 + b * r33;
  b21 = r21;           b22 = r22;           b23 = r23;
  b31 = -b * r11 + a * r31; b32 = -b * r12 + a * r32; b33 = -b * r13 + a * r33;

  // third givens rotation (ch,sh,0,0)
  QRGivensQuaternion(b22, b32, ch3, sh3);
  a = 1 - 2 * sh3 * sh3;
  b = 2 * ch3 * sh3;
  // R is now set to desired value
  r11 = b11;             r12 = b12;           r13 = b13;
  r21 = a * b21 + b * b31;     r22 = a * b22 + b * b32;   r23 = a * b23 + b * b33;
  r31 = -b * b21 + a * b31;    r32 = -b * b22 + a * b32;  r33 = -b * b23 + a * b33;

  // construct the cumulative rotation Q=Q1 * Q2 * Q3
  // the number of floating point operations for three quaternion multiplications
  // is more or less comparable to the explicit form of the joined matrix.
  // certainly more memory-efficient!
  float sh12 = sh1 * sh1;
  float sh22 = sh2 * sh2;
  float sh32 = sh3 * sh3;

  q11 = (-1 + 2 * sh12) * (-1 + 2 * sh22);
  q12 = 4 * ch2 * ch3 * (-1 + 2 * sh12) * sh2 * sh3 + 2 * ch1 * sh1 * (-1 + 2 * sh32);
  q13 = 4 * ch1 * ch3 * sh1 * sh3 - 2 * ch2 * (-1 + 2 * sh12) * sh2 * (-1 + 2 * sh32);

  q21 = 2 * ch1 * sh1 * (1 - 2 * sh22);
  q22 = -8 * ch1 * ch2 * ch3 * sh1 * sh2 * sh3 + (-1 + 2 * sh12) * (-1 + 2 * sh32);
  q23 = -2 * ch3 * sh3 + 4 * sh1 * (ch3 * sh1 * sh3 + ch1 * ch2 * sh2 * (-1 + 2 * sh32));

  q31 = 2 * ch2 * sh2;
  q32 = 2 * ch3 * (1 - 2 * sh22) * sh3;
  q33 = (-1 + 2 * sh22) * (-1 + 2 * sh32);
}

inline
void svd(// input A
  float* a,
  //		float a11, float a12, float a13,
  //		float a21, float a22, float a23,
  //		float a31, float a32, float a33,
  // output U
  float* u,
  //		float &u11, float &u12, float &u13,
  //		float &u21, float &u22, float &u23,
  //		float &u31, float &u32, float &u33,
  // output S
  float* s,
  //		float &s11, float &s12, float &s13,
  //		float &s21, float &s22, float &s23,
  //        float &s31, float &s32, float &s33,
  // output V
  float* v
  //		float &v11, float &v12, float &v13,
  //		float &v21, float &v22, float &v23,
  //		float &v31, float &v32, float &v33
) {
  // normal equations matrix
  //  float ATA11, ATA12, ATA13;
  //  float ATA21, ATA22, ATA23;
  //  float ATA31, ATA32, ATA33;
  float ATA[9];

  multAtB(a, a, ATA);

  // symmetric eigenalysis
  float qV[4];
  jacobiEigenanlysis(
    ATA[0], ATA[3],  ATA[4],  ATA[6],  ATA[7],  ATA[8], qV);
  quatToMat3(qV, v);//v11, v12, v13, v21, v22, v23, v31, v32, v33);

  //  float b11, b12, b13;
  //  float b21, b22, b23;
  //  float b31, b32, b33;
  float b[9];
  multAB(a, v, b);

  // sort singular values and find V
  sortSingularValues(b[0], b[1], b[2],
                     b[3], b[4], b[5],
                     b[6], b[7], b[8],

                     v[0], v[1], v[2],
                     v[3], v[4], v[5],
                     v[6], v[7], v[8]);

  // QR decomposition
  QRDecomposition(b[0], b[1], b[2],
                  b[3], b[4], b[5],
                  b[6], b[7], b[8],

                  u[0], u[1], u[2],
                  u[3], u[4], u[5],
                  u[6], u[7], u[8],

                  s[0], s[1], s[2],
                  s[3], s[4], s[5],
                  s[6], s[7], s[8]);
}

/// polar decomposition can be reconstructed trivially from SVD result
// A = UP
//void pd(float a11, float a12, float a13,
//        float a21, float a22, float a23,
//        float a31, float a32, float a33,
//        // output U
//        float &u11, float &u12, float &u13,
//        float &u21, float &u22, float &u23,
//        float &u31, float &u32, float &u33,
//        // output P
//        float &p11, float &p12, float &p13,
//        float &p21, float &p22, float &p23,
//        float &p31, float &p32, float &p33)
//{
//    float w11, w12, w13, w21, w22, w23, w31, w32, w33;
//    float s11, s12, s13, s21, s22, s23, s31, s32, s33;
//    float v11, v12, v13, v21, v22, v23, v31, v32, v33;
//
//    svd(a11, a12, a13, a21, a22, a23, a31, a32, a33,
//        w11, w12, w13, w21, w22, w23, w31, w32, w33,
//        s11, s12, s13, s21, s22, s23, s31, s32, s33,
//        v11, v12, v13, v21, v22, v23, v31, v32, v33);
//
//    // P = VSV'
//    float t11, t12, t13, t21, t22, t23, t31, t32, t33;
//    multAB(v11, v12, v13, v21, v22, v23, v31, v32, v33,
//           s11, s12, s13, s21, s22, s23, s31, s32, s33,
//           t11, t12, t13, t21, t22, t23, t31, t32, t33);
//
//    multAB(t11, t12, t13, t21, t22, t23, t31, t32, t33,
//           v11, v21, v31, v12, v22, v32, v13, v23, v33,
//           p11, p12, p13, p21, p22, p23, p31, p32, p33);
//
//    // U = WV'
//    multAB(w11, w12, w13, w21, w22, w23, w31, w32, w33,
//           v11, v21, v31, v12, v22, v32, v13, v23, v33,
//           u11, u12, u13, u21, u22, u23, u31, u32, u33);
//}

}
namespace nr {
const float kEpsilon = std::numeric_limits<float>::epsilon();
//const float kEpsilon = 1e-20f;

inline float pythag(float a, float b) {
  float absa, absb;
  absa = fabs(a);
  absb = fabs(b);
  if (absa > absb) return absa * sqrtf(1.0f + dj::Square(absb / absa));
  else return (absb == 0.0f ? 0.0f : absb * sqrtf(1.0f + dj::Square(absa / absb)));
}

template <int m, int n>
inline bool svdcmp(float (*u)[n], float w[], float (*v)[n]) {
  bool flag;
  int i, its, j, jj, k, l, nm;
  float anorm, c, f, g, h, s, scale, x, y, z;
  float rv1[n];
  g = scale = anorm = 0.0f;
  bool ret = true;
  for (i = 0; i < n; i++) {
    l = i + 2;
    rv1[i] = scale * g;
    g = s = scale = 0.0f;
    if (i < m) {
      for (k = i; k < m; k++) scale += dj::Abs(u[k][i]);
      if (scale != 0.0) {
        for (k = i; k < m; k++) {
          u[k][i] /= scale;
          s += u[k][i] * u[k][i];
        }
        f = u[i][i];
        g = -NR_SIGN(sqrt(s), f);
        h = f * g - s;
        u[i][i] = f - g;
        for (j = l - 1; j < n; j++) {
          for (s = 0.0f, k = i; k < m; k++) s += u[k][i] * u[k][j];
          f = s / h;
          for (k = i; k < m; k++) u[k][j] += f * u[k][i];
        }
        for (k = i; k < m; k++) u[k][i] *= scale;
      }
    }
    w[i] = scale * g;
//    printf("%d %f %f %f\n", i, w[i], scale, g);
    g = s = scale = 0.0f;
    if (i + 1 <= m && i + 1 != n) {
      for (k = l - 1; k < n; k++) scale += dj::Abs(u[i][k]);
      if (scale != 0.0) {
        for (k = l - 1; k < n; k++) {
          u[i][k] /= scale;
          s += u[i][k] * u[i][k];
        }
        f = u[i][l - 1];
        g = -NR_SIGN(sqrt(s), f);
        h = f * g - s;
        u[i][l - 1] = f - g;
        for (k = l - 1; k < n; k++) rv1[k] = u[i][k] / h;
        for (j = l - 1; j < m; j++) {
          for (s = 0.0f, k = l - 1; k < n; k++) s += u[j][k] * u[i][k];
          for (k = l - 1; k < n; k++) u[j][k] += s * rv1[k];
        }
        for (k = l - 1; k < n; k++) u[i][k] *= scale;
      }
    }
    anorm = dj::Max(anorm, (dj::Abs(w[i]) + dj::Abs(rv1[i])));
  }
  for (i = n - 1; i >= 0; i--) {
    if (i < n - 1) {
      if (g != 0.0f) {
        for (j = l; j < n; j++)
          v[j][i] = (u[i][j] / u[i][l]) / g;
        for (j = l; j < n; j++) {
          for (s = 0.0f, k = l; k < n; k++) s += u[i][k] * v[k][j];
          for (k = l; k < n; k++) v[k][j] += s * v[k][i];
        }
      }
      for (j = l; j < n; j++) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0f;
    g = rv1[i];
    l = i;
  }
  for (i = dj::Min(m, n) - 1; i >= 0; i--) {
    l = i + 1;
    g = w[i];
    for (j = l; j < n; j++) u[i][j] = 0.0f;
    if (g != 0.0f) {
      g = 1.0f / g;
      for (j = l; j < n; j++) {
        for (s = 0.0f, k = l; k < m; k++) s += u[k][i] * u[k][j];
        f = (s / u[i][i]) * g;
        for (k = i; k < m; k++) u[k][j] += f * u[k][i];
      }
      for (j = i; j < m; j++) u[j][i] *= g;
    } else for (j = i; j < m; j++) u[j][i] = 0.0f;
    ++u[i][i];
  }
  for (k = n - 1; k >= 0; k--) {
    for (its = 0; its < 30; its++) {
      flag = true;
      for (l = k; l >= 0; l--) {
        nm = l - 1;
        if (l == 0 || dj::Abs(rv1[l]) <= kEpsilon * anorm) {
          flag = false;
          break;
        }
        if (dj::Abs(w[nm]) <= kEpsilon * anorm) break;
      }
      if (flag) {
        c = 0.0f;
        s = 1.0f;
        for (i = l; i < k + 1; i++) {
          f = s * rv1[i];
          rv1[i] = c * rv1[i];
          if (dj::Abs(f) <= kEpsilon * anorm) break;
          g = w[i];
          h = pythag(f, g);
          w[i] = h;
          h = 1.0f / h;
          c = g * h;
          s = -f * h;
          for (j = 0; j < m; j++) {
            y = u[j][nm];
            z = u[j][i];
            u[j][nm] = y * c + z * s;
            u[j][i] = z * c - y * s;
          }
        }
      }
      z = w[k];
      if (l == k) {
        if (z < 0.0f) {
          w[k] = -z;
          for (j = 0; j < n; j++) v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == 29) {
        ret = false;
        std::cerr << ("no convergence in 30 svdcmp iterations") << std::endl;
      }
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0f * h * y);
      g = pythag(f, 1.0f);
      f = ((x - z) * (x + z) + h * ((y / (f + NR_SIGN(g, f))) - h)) / x;
      c = s = 1.0;
      for (j = l; j <= nm; j++) {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = pythag(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y *= c;
        for (jj = 0; jj < n; jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x * c + z * s;
          v[jj][i] = z * c - x * s;
        }
        z = pythag(f, h);
        w[j] = z;
        if (z) {
          z = 1.0f / z;
          c = f * z;
          s = h * z;
        }
        f = c * g + s * y;
        x = c * y - s * g;
        for (jj = 0; jj < m; jj++) {
          y = u[jj][j];
          z = u[jj][i];
          u[jj][j] = y * c + z * s;
          u[jj][i] = z * c - y * s;
        }
      }
      rv1[l] = 0.0f;
      rv1[k] = f;
      w[k] = x;
    }
  }
  return ret;
}
}
#undef NRANSI
#endif // SVD_H
