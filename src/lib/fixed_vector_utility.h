#ifndef FIXED_VECTOR_UTILITY_H_
#define FIXED_VECTOR_UTILITY_H_
#include "utility_function.h"
#include "macro_constant.h"
#include "print_macro.h"
//#define NOMINMAX
#undef max
#undef min
namespace dj {

template <class T>
inline void SubVec2(const T* a, const T* b, T* result) {
  result[0] = a[0] - b[0];
  result[1] = a[1] - b[1];
}

template <class T>
inline void SubVec3(const T* a, const T* b, T* result) {
  result[0] = a[0] - b[0];
  result[1] = a[1] - b[1];
  result[2] = a[2] - b[2];
}

template <class T>
inline void AddVec2(const T* a, const T* b, T* result) {
  result[0] = a[0] + b[0];
  result[1] = a[1] + b[1];
}

template <class T>
inline void AddVec3(const T* a, const T* b, T* result) {
  result[0] = a[0] + b[0];
  result[1] = a[1] + b[1];
  result[2] = a[2] + b[2];
}

template <class T>
inline T Normalize2(T* vector) {
  T norm =  sqrt(vector[0] * vector[0] + vector[1] * vector[1]) ;
  if (norm < EPSILON) {
    return 0;
  }
  T invNorm = 1 / norm ;
  vector[0] *= invNorm ;
  vector[1] *= invNorm ;
  return norm;
}

template <class T>
inline T Normalize2(T& x, T& y) {
  T norm = hypot(x, y);
  if (norm < EPSILON) {
    return 0 ;
  }
  T invNorm = 1 / norm ;
  x *= invNorm ;
  y *= invNorm ;
  return norm ;
}

template <class T>
inline T Normalize3(T* a) {
  T norm =  sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) ;
  //  if (norm < EPSILON) {
  //    return 0 ;
  //  }
  T invNorm = 1 / norm ;
  if (invNorm != invNorm) {
    std::cerr << CURRENT_LINE << "=> NaN when normalizing!!" << std::endl;
    return 0;
  }
  a[0] *= invNorm;
  a[1] *= invNorm;
  a[2] *= invNorm;
  return norm;
}

template <class T>
inline T Norm2(T* v) {
  return sqrt(v[0] * v[0] + v[1] * v[1]);
} // #Norm2#

template <class T>
inline T Norm3(T* v) {
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
} // #Norm3#


/// Compute v1 X v2
template <class T>
inline void Cross3(const T* v1, const T* v2, T* cross) {
  cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
  cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
  cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

template <class T>
inline T Cross2(const T* v1, const T* v2) {
  return v1[0] * v2[1] - v1[1] * v2[0];
}

template <class T>
inline T Dot2(const T* v1, const T* v2) {
  return v1[0] * v2[0] + v1[1] * v2[1];
}

template <class T>
inline T Dot3(const T* v1, const T* v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

template <class T>
inline T Distance2(T* a, T *b) {
  return hypot(a[0] - b[0], a[1] - b[1]) ;
}

template <class T>
inline T Distance3(const T* a, const T* b) {
  return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                   (a[2] - b[2]) * (a[2] - b[2]));
}

/*
 From http://www.dennis2society.de/main/painless-tetrahedral-barycentric-mapping
  Calculate the determinant for a 4x4 matrix based on this example:
  http://www.euclideanspace.com/maths/algebra/matrix/functions/determinant/fourD/index.htm
  This function takes four Vec4f as row vectors and calculates the resulting matrix' determinant
  using the Laplace expansion.
const float Determinant4x4( const Vec4f& v0,
                            const Vec4f& v1,
                            const Vec4f& v2,
                            const Vec4f& v3 )
{
float det = v0.w*v1.z*v2.y*v3.x - v0.z*v1.w*v2.y*v3.x -
            v0.w*v1.y*v2.z*v3.x + v0.y*v1.w*v2.z*v3.x +

            v0.z*v1.y*v2.w*v3.x - v0.y*v1.z*v2.w*v3.x -
            v0.w*v1.z*v2.x*v3.y + v0.z*v1.w*v2.x*v3.y +

            v0.w*v1.x*v2.z*v3.y - v0.x*v1.w*v2.z*v3.y -
            v0.z*v1.x*v2.w*v3.y + v0.x*v1.z*v2.w*v3.y +

            v0.w*v1.y*v2.x*v3.z - v0.y*v1.w*v2.x*v3.z -
            v0.w*v1.x*v2.y*v3.z + v0.x*v1.w*v2.y*v3.z +

            v0.y*v1.x*v2.w*v3.z - v0.x*v1.y*v2.w*v3.z -
            v0.z*v1.y*v2.x*v3.w + v0.y*v1.z*v2.x*v3.w +

            v0.z*v1.x*v2.y*v3.w - v0.x*v1.z*v2.y*v3.w -
            v0.y*v1.x*v2.z*v3.w + v0.x*v1.y*v2.z*v3.w;
    return det;
}
*/
template <class T>
inline T Determinant4(const T* p0, const T* p1, const T* p2, const T* p3) {
  T det = p1[2] * p2[1] * p3[0] - p0[2] * p2[1] * p3[0] -
          p1[1] * p2[2] * p3[0] + p0[1] * p2[2] * p3[0] +

          p0[2] * p1[1] * p3[0] - p0[1] * p1[2] * p3[0] -
          p1[2] * p2[0] * p3[1] + p0[2] * p2[0] * p3[1] +

          p1[0] * p2[2] * p3[1] - p0[0] * p2[2] * p3[1] -
          p0[2] * p1[0] * p3[1] + p0[0] * p1[2] * p3[1] +

          p1[1] * p2[0] * p3[2] - p0[1] * p2[0] * p3[2] -
          p1[0] * p2[1] * p3[2] + p0[0] * p2[1] * p3[2] +

          p0[1] * p1[0] * p3[2] - p0[0] * p1[1] * p3[2] -
          p0[2] * p1[1] * p2[0] + p0[1] * p1[2] * p2[0] +

          p0[2] * p1[0] * p2[1] - p0[0] * p1[2] * p2[1] -
          p0[1] * p1[0] * p2[2] + p0[0] * p1[1] * p2[2];
  return det;
}

template <class T>
inline T Determinant3(const T (*m)[3]) {
  return + m[0][0] * m[1][1] * m[2][2]
         + m[0][1] * m[1][2] * m[2][0]
         + m[0][2] * m[1][0] * m[2][1]
         - m[0][0] * m[1][2] * m[2][1]
         - m[0][1] * m[1][0] * m[2][2]
         - m[0][2] * m[1][1] * m[2][0];
}

template <class T>
inline void Inverse3(const T (*m)[3], T (*inv)[3]) {
  P(Determinant3(m));
  T invdet = 1 / Determinant3(m);
  P(invdet);
  inv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
  inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
  inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
  inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
  inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
  inv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
  inv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
  inv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
  inv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
}

/**
 * @param p_i: the tetrahedron vertex position
 * @param p: the point barycentric coordinate of which this function computes
 * @param barycentricCoord: result of barycentric coordinate computation
 * @return distance from p to tet's mass center if point is out side of tet
 */
template <class T>
inline T GetBarycentricCoordinate(const T* p0, const T* p1, const T* p2, const T* p3,
                                  const T* p, T* barycentricCoord) {
  const T det = Determinant4(p0, p1, p2, p3);
  if (det < EPSILON && det > -EPSILON) {
    std::cerr << "GetBarycentricCoordinate() => Tet is almost flat!!" << std::endl;
    return T(1e10);
  }
  const T det_inv = T(1.0) / det;
  barycentricCoord[0] = Determinant4(p, p1, p2, p3) * det_inv;
  barycentricCoord[1] = Determinant4(p0, p, p2, p3) * det_inv;
  barycentricCoord[2] = Determinant4(p0, p1, p, p3) * det_inv;
  barycentricCoord[3] = Determinant4(p0, p1, p2, p) * det_inv;
  if (barycentricCoord[0] < -EPSILON || barycentricCoord[1] < -EPSILON
      || barycentricCoord[2] < -EPSILON || barycentricCoord[3] < -EPSILON) {
    T center[3] = {(p0[0] + p1[0] + p2[0] + p3[0]) * T(0.25),
                   (p0[1] + p1[1] + p2[1] + p3[1]) * T(0.25),
                   (p0[2] + p1[2] + p2[2] + p3[2]) * T(0.25),
                  };
    return Distance3(center, p);
  } else {
    return T(-1.0);
  }
}

/**
 * @param p_i: the tetrahedron vertex position
 * @param p: the point barycentric coordinate of which this function computes
 * @param barycentricCoord: result of barycentric coordinate computation
 * @return  true if p is inside the tet, false o/w
 */
//#include "print_macro.h"
template <class T>
inline bool IsInsideTet(const T* p0, const T* p1, const T* p2, const T* p3,
                        const T* p, float threshold, T* barycentricCoord) {
  const float kEpsilon = 1e-12f;
  const T det = Determinant4(p0, p1, p2, p3);
  //  if (det < EPSILON && det > -EPSILON) {
  if (det < kEpsilon && det > -kEpsilon) {
    std::cerr << "GetBarycentricCoordinate() => Tet volume (" << det << ") is almost flat!!" << std::endl;
    return false;
  }
  const T det_inv = T(1.0) / det;
  barycentricCoord[0] = Determinant4(p, p1, p2, p3) * det_inv;
  barycentricCoord[1] = Determinant4(p0, p, p2, p3) * det_inv;
  barycentricCoord[2] = Determinant4(p0, p1, p, p3) * det_inv;
  barycentricCoord[3] = Determinant4(p0, p1, p2, p) * det_inv;

  //  if (barycentricCoord[0] < -kEpsilon
  //      || barycentricCoord[1] < -kEpsilon
  //      || barycentricCoord[2] < -kEpsilon
  //      || barycentricCoord[3] < -kEpsilon) {
  if (barycentricCoord[0] < -threshold
      || barycentricCoord[1] < -threshold
      || barycentricCoord[2] < -threshold
      || barycentricCoord[3] < -threshold) {
    return false;
  } else {
    return true;
  }
}


template <class T>
inline T GetBarycentricCoordinate(const T* p0, const T* p1, const T* p2, const T* p3,
                                  const T det, const T* center,
                                  const T* p, T* barycentricCoord) {
  if (det < EPSILON && det > -EPSILON) {
    return T(1e10);
  }
  const T det_inv = T(1.0) / det;
  barycentricCoord[0] = Determinant4(p, p1, p2, p3) * det_inv;
  barycentricCoord[1] = Determinant4(p0, p, p2, p3) * det_inv;
  barycentricCoord[2] = Determinant4(p0, p1, p, p3) * det_inv;
  barycentricCoord[3] = Determinant4(p0, p1, p2, p) * det_inv;
  if (barycentricCoord[0] < -EPSILON || barycentricCoord[1] < -EPSILON
      || barycentricCoord[2] < -EPSILON || barycentricCoord[3] < -EPSILON) {
    return Distance3(p, center);
  } else {
    return T(-1.0);
  }
}


/// compute barycentric coordinate of p for triangle (p0, p1, p2)
template <class T>
inline void GetBarycentricCoordinate(const T* p0, const T* p1, const T* p2,
                                     const T* p, T* barycentricCoord) {
  T vec[3][3] = {
    p[0] - p0[0], p[1] - p0[1], p[2] - p0[2],
    p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2],
    p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2],
  };
  float fake_normal[3];
  dj::Cross3(vec[1], vec[2], fake_normal);
  float denominator = dj::Dot3(fake_normal, fake_normal);
  denominator = T(1.0) / denominator;

  float nominator[2][3];
  dj::Cross3(vec[1], vec[0], nominator[0]);
  dj::Cross3(vec[0], vec[2], nominator[1]);

  barycentricCoord[1] = dj::Dot3(fake_normal, nominator[1]) * denominator;
  barycentricCoord[2] = dj::Dot3(fake_normal, nominator[0]) * denominator;
  barycentricCoord[0] = 1 - barycentricCoord[1] - barycentricCoord[2];
}

template <typename T1, typename T2>
struct MulType {
private:
  static T1* a;
  static T2* b;
public:
  typedef decltype((*a) * (*b)) Type;
};

template <int d, class T>
struct FixedVecOp {
  inline static T Norm2Square(const T* x) {
    T sum = (T) 0;
    for (int i = 0; i < d; ++i) {
      sum += (x[i] * x[i]);
    }
    return sum;
  }

  inline static T MagnitudeSquare(const T* x) {
    T sum = T(0);
    for (int i = 0; i < d; ++i) {
      sum += (x[i] * x[i]);
    }
    return sum;
  }

  inline static T Magnitude(const T* x) {
    T sum = T(0);
    for (int i = 0; i < d; ++i) {
      sum += (x[i] * x[i]);
    }
    return std::sqrt(sum);
  }

  inline static T Normalize(T* x) {
    T length = FixedVecOp<d, T>::Magnitude(x);
    if (length < EPSILON) {
      std::cerr << "FixedVectorOp::Normalize() => try to normalize 0 vector." << std::endl;
      int* a = NULL;
      *a = 0;
      assert(false);
      return T(0);
    } else {
      T length_inv = T(1.0) / length;
      for (int i = 0; i < d; ++i) {
        x[i] *= length_inv;
      }
      return length;
    }
  }

  inline static bool Equal(const T* x1, const T* x2) {
    for (int i = 0; i < d; ++i) {
      if (x1[i] != x2[i]) return false;
    }
    return true;
  }

  inline static bool LessThan(const T* x1, const T* x2) {
    for (int i = 0; i < d; ++i) {
      if (x1[i] > x2[i]) return false;
      else if (x1[i] < x2[i]) return true;
    }
    return false;
  }

  inline static void Copy(const T* from, T* to) {
    for (int i = 0; i < d; ++i) {
      to[i] = from[i];
    }
  }

  inline static void Scale(T* x, T scale) {
    for (int i = 0; i < d; ++i) {
      x[i] *= scale;
      //      x[i] = scale * x[i];
    }
  }

  inline static void Minus(const T* x, const T* y, T* result) {
    for (int i = 0; i < d; ++i) {
      result[i] = x[i] - y[i];
    }
  }

  inline static void Plus(const T* x, const T* y, T* const result) {
    for (int i = 0; i < d; ++i) {
      result[i] = x[i] + y[i];
    }
  }

  template <class T1>
  inline static typename MulType<T, T1>::Type Dot(const T* x, const T1* y) {
    typedef typename MulType<T, T1>::Type ResultType;
    ResultType dot = ResultType();
    for (int i = 0; i < d; ++i) {
      dot += x[i] * y[i];
    }
    return dot;
  }

  inline static T DistanceSquare(const T* x, const T* y) {
    T sum = T(0);
    for (int i = 0; i < d; ++i) {
      sum += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sum;
  }

  inline static T Distance(const T* x, const T* y) {
    T sum = T(0);
    for (int i = 0; i < d; ++i) {
      sum += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return std::sqrt(sum);
  }

  inline static void Fill(T* x, const T value) {
    for (int i = 0; i < d; ++i) {
      x[i] = value;
    }
  }

  inline static T MinAbs(const T* x) {
    T min = T(1e8);
    for (int i = 0; i < d; i++) {
      T abs = dj::Abs(x[i]);
      if (min > abs) {
        min = abs;
      }
    }
    return min;
  }

  inline static T MaxAbs(const T* x) {
    T max = T(EPSILON);
    for (int i = 0; i < d; i++) {
      T abs = dj::Abs(x[i]);
      if (max < abs) {
        max = abs;
      }
    }
    return max;
  }

  inline static T Min(const T* x) {
    T min = T(1e8);
    for (int i = 0; i < d; i++) {
      if (min > x[i]) {
        min = x[i];
      }
    }
    return min;
  }

  inline static T Max(const T* x) {
    T max = std::numeric_limits<T>::min();
    for (int i = 0; i < d; i++) {
      if (max < x[i]) {
        max = x[i];
      }
    }
    return max;
  }
};

}
#endif // FIXED_VECTOR_UTILITY_H_
