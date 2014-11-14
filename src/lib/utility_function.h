#ifndef UTILITY_FUNCTION_H_
#define UTILITY_FUNCTION_H_
#include <type_traits>
#include <cmath>
#include "macro_constant.h"
namespace dj
{
template <class T>
inline void Swap(T& a, T& b)
{
  T temp = a;
  a = b ;
  b = temp ;
}
inline bool Eq(float a, float b, float epsilon = EPSILON)
{
  float diff = a - b;
  return (diff > -epsilon && diff < epsilon);
}

inline bool Eq(double a, double b, double epsilon = EPSILON)
{
  double diff = a - b;
  return (diff > -epsilon && diff < epsilon);
}


inline float Degree2Radian(float degree) {
  return degree / 180.0f * PI;
}

inline double Degree2Radian(double degree) {
  return degree / 180.0 * PI;
}

inline float Radian2Degree(float radian)
{
  return radian / PI * 180.0f;
}

inline double Radian2Degree(double radian)
{
  return radian / PI * 180.0;
}

template <class T>
inline T Square(T x)
{
  return x * x;
}

template <class T>
inline T Abs(T x)
{
  return (x >= 0) ? x : -x;
}

template <class T>
inline T Max(T a, T b)
{
  return (a > b) ? a : b;
}

template <class T>
inline T Min(T a, T b)
{
  return (a <= b) ? a : b;
}

template <class T>
inline T Sign(T x)
{
  return (x >= T(0)) ? T(1) : T(-1);
}

template <class T>
inline bool NaN(T x)
{
  return x == x;
}

//------------------------------------------------------------------------------
// Interpolation functions
//------------------------------------------------------------------------------
template <class T>
inline T BilinearInterpolate(const T lower_left,
                             const T lower_right,
                             const T upper_left,
                             const T upper_right,
                             const float distance2left,
                             const float distance2lower)
{
  return distance2left * (distance2lower * upper_right + (1 - distance2lower) * lower_right) +
         (1 - distance2left) * (distance2lower * upper_left + (1 - distance2lower) * lower_left);
}


template <class T>
inline T Tex2D(const T* texture,
               int nx, int ny,
               T x, T y)
{
  (void) ny;
  int x0 = int(x);
  int y0 = int(y);
  int lower_left = x0 + y0 * nx;
  return BilinearInterpolate(texture[lower_left],
                             texture[lower_left + 1],
                             texture[lower_left + nx],
                             texture[lower_left + nx + 1],
                             x - x0, y - y0);
}

template <class T>
inline T tex2DNormalized(const T* texture,
                         int nx, int ny,
                         T x, T y)
{
  return Tex2D(texture, nx, ny, x * nx, y * ny);
}

template <typename T, typename T2>
inline T LinearInterpolate(const T left, const T right, T2 distance2left)
{
  return distance2left * right + (1 - distance2left) * left;
}

template <class T>
inline T Tex1D(const T* texture,
               int nx,
               T x)
{
  (void) nx;
  int x0 = int(x);
  return LinearInterpolate(texture[x0], texture[x0 + 1], x - x0);
}

template <class T>
inline T Tex1DNormalized(T* texture,
                         int nx,
                         T x)
{
  return Tex1D(texture, nx, x * nx);
}



} // namespace
#endif

