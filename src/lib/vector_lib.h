#ifndef VECTOR_LIB_H_
#define VECTOR_LIB_H_
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <sstream>
#include <cstring>
#include <vector>
#include <limits>
#include "fixed_vector_utility.h"
#include "fixed_vector.h"
#include "utility_function.h"
#include "fixed_matrix_utility.h"
#include "fixed_matrix.h"
#include "macro_constant.h"
//#include "print_macro.h"
#define NaN(x) ((x) != (x))

namespace dj
{
template <class T>
std::string Vec3ToString(T* vec)
{
  std::stringstream str;
  str << "[" << vec[0] << ", " << vec[1] << ", " << vec[2] << "]";
  return str.str();
}

template <class T>
inline std::vector<T> make_vector(T* vec, int size)
{
  return std::vector<T>(vec, vec + size);
}

template <class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec)
{
  if (vec.size() == 0) {
    out << "{ }";
  } else {
    out << "{" << vec[0];
    for (unsigned int i = 1; i < vec.size(); i++) {
      out << ", " << vec[i];
    }
    out << "}";
  }
  return out;
}

//------------------------------------------------------------------------------
// Matrix library


template <class T, int row_num, int col_num>
std::string Array2String(const T (*matrix)[col_num])
{
  if (row_num < 1 || col_num < 1) {
    return "[]";
  } else {
    std::stringstream out;
    out << "[";
    int i = 0;
    for (; i < row_num - 1; ++i) {
      for (int j = 0; j < col_num - 1; j++) {
        out << matrix[i][j] << " ";
      }
      out << matrix[i][col_num - 1] << "; ";
    }
    i = row_num - 1;
    for (int j = 0; j < col_num - 1; j++) {
      out << matrix[i][j] << " ";
    }
    out << matrix[i][col_num - 1];
    out << "]";
    return out.str();
  }
}

template <class T>
T ManhattanDistance2(T* a, T* b)
{
  return dj::Abs(a[0] - b[0]) + dj::Abs(a[1] - b[1]);
}

template <class T>
T ManhattanDistance3(T* a, T* b)
{
  return dj::Abs(a[0] - b[0]) + dj::Abs(a[1] - b[1]) + dj::Abs(a[2] - b[2]);
}

template <class T, int row_num, int col_num>
std::ostream& operator<<(std::ostream& out, const FixedMatrixWrapper<T, row_num, col_num>& matrix)
{
  out << Array2String<T, row_num, col_num>(matrix.matrix_);
  return out;
}

template <class T, int row_num, int col_num>
std::ostream& operator<<(std::ostream& out, const FixedMatrix<T, row_num, col_num>& matrix)
{
  out << Array2String<T, row_num, col_num>((T (*)[col_num]) matrix.matrix_);
  return out;
}

template <class T, int row_num, int col_num>
inline const std::pair<int, int>& Size(const FixedMatrixWrapper<T, row_num, col_num>&)
{
  static const std::pair<int, int> size(row_num, col_num);
  return size;
}

template <class T, int row_num, int col_num>
inline const std::pair<int, int>& Size(const FixedMatrix<T, row_num, col_num>&)
{
  static const std::pair<int, int> size(row_num, col_num);
  return size;
}

template <class T, int d, bool wrapped>
inline int Size(const Vec<T, d, wrapped>&)
{
  return Vec<T, d, wrapped>::kDimension;
}



// End of matrix library
//------------------------------------------------------------------------------
template <class T>
bool SegmentTriangleIntersection(T *tri[3], T* p0, T* p1)
{
  using namespace dj;
  T vec[2][3];
  SubVec3(tri[1], tri[0], vec[0]);
  dj::SubVec3(tri[2], tri[0], vec[1]);
  T segment_direction[3];
  dj::SubVec3(p1, p0, segment_direction);
  T tri_normal[3];
  dj::Cross3(vec[0], vec[1], tri_normal);
  dj::Normalize3(tri_normal);
  T v0_to_p0[3];
  dj::SubVec3(p0, tri[0], v0_to_p0);
  T dots[2] = {
    dj::Dot3(tri_normal, v0_to_p0),
    dj::Dot3(tri_normal, segment_direction),
  };
  if (dots[0] * dots[1] > 0) return false;
  T v1_to_v2[3];
  dj::SubVec3(tri[2], tri[1], v1_to_v2);

  auto InSideTriangle = [&](T * p) -> bool {
    T v0_to_p[3];
    dj::SubVec3(p, tri[0], v0_to_p);
    T edge[3];
    dj::Cross3(vec[0], v0_to_p, edge);
    if (dj::Dot3(edge, tri_normal) < 0) return false;
    dj::Cross3(v0_to_p, vec[1], edge);
    if (dj::Dot3(edge, tri_normal) < 0) return false;

    T v1_to_p[3];
    dj::SubVec3(p, tri[1], v1_to_p);
    dj::Cross3(v1_to_v2, v1_to_p, edge);
    if (dj::Dot3(edge, tri_normal) < 0) return false;
    return true;
  };

  if (dj::Abs(dots[1]) < 1e-6f) {
    // coplanar
    return (dj::Abs(dots[0]) < 1e-6f) // distance to plance is 0
           && (InSideTriangle(p0) || InSideTriangle(p1)); // one point inside and one outside
  } else {
    T intersection[3];
    if (dj::Abs(dots[1]) < 1e-6f) {
      intersection[0] = p0[0];
      intersection[1] = p0[1];
      intersection[2] = p0[2];
    } else {
      T ratio = dots[0] / dots[1];
      if (ratio < -1) return false; // out side of segment
      intersection[0] = p0[0] - ratio * segment_direction[0];
      intersection[1] = p0[1] - ratio * segment_direction[1];
      intersection[2] = p0[2] - ratio * segment_direction[2];
    }
    bool intersect = InSideTriangle(intersection);
    if (intersect) {
      P(Vec3ToString(intersection));
      P(Vec3ToString(tri[0]));
      P(Vec3ToString(tri[1]));
      P(Vec3ToString(tri[2]));
      P(Vec3ToString(p0));
      P(Vec3ToString(p1));
    }
    return intersect;
  }
}


template <class T>
inline T Point2SegmentDistance(const T* p0, const T* p1, const T* p, T& projection_parameter)
{
  T vec1[3];
  T vec2[3];
  dj::SubVec3(p1, p0, vec1);
  T length = Normalize3(vec1);
  dj::SubVec3(p, p0, vec2);
  T projection_length = Dot3(vec1, vec2);
  T distance = 0;
  if (projection_length < 0) {
    distance = Distance3(p0, p);
    projection_parameter = 0;
  } else if (projection_length > length) {
    distance = Distance3(p1, p);
    projection_parameter = 1;
  } else {
    float projection_point[3] = {
      vec1[0] * projection_length + p0[0],
      vec1[1] * projection_length + p0[1],
      vec1[2] * projection_length + p0[2]
    };
    distance = Distance3(projection_point, p);
    projection_parameter = projection_length / length;
  }
  return distance;
}

template <class T>
inline T PointLineDistance(const T* p0, const T* p1, const T* p, T* length = NULL)
{
  T vec1[3];
  T vec2[3];
  dj::SubVec3(p1, p0, vec1);
  Normalize3(vec1);
  dj::SubVec3(p, p0, vec2);
  T projection_length = Dot3(vec1, vec2);
  if (length) {
    *length = projection_length;
  }
  T distance = 0;
  distance = std::sqrt(dj::Square(Norm3(vec2)) - dj::Square(projection_length));
  return distance;
}


template <class T>
inline void Translate(T* position_array, int num, T dx, T dy, T dz)
{
  for (int i = 0, v3 = 0; i < num; ++i, v3 += 3) {
    position_array[v3 + 0] += dx;
    position_array[v3 + 1] += dy;
    position_array[v3 + 2] += dz;
  }
}

template <class T>
inline void GetCenter(const T* position_array, int num, T* center)
{
  center[0] = center[1] = center[2] = 0;
  for (int i = 0, v3 = 0; i < num; ++i, v3 += 3) {
    center[0] += position_array[v3 + 0];
    center[1] += position_array[v3 + 1];
    center[2] += position_array[v3 + 2];
  }
  center[0] /= num;
  center[1] /= num;
  center[2] /= num;
}

template <class T>
inline void Scale(T* position_array, int num, T sx, T sy, T sz)
{
  for (int i = 0, v3 = 0; i < num; ++i, v3 += 3) {
    position_array[v3 + 0] *= sx;
    position_array[v3 + 1] *= sy;
    position_array[v3 + 2] *= sz;
  }
}

template <class T>
inline void ScaleAroundCenter(T* position_array, int num, T sx, T sy, T sz)
{
  float center[3];
  dj::GetCenter(position_array, num, center);
  for (int i = 0, v3 = 0; i < num; ++i, v3 += 3) {
    position_array[v3 + 0] = center[0] + (position_array[v3 + 0] - center[0]) * sx;
    position_array[v3 + 1] = center[1] + (position_array[v3 + 1] - center[1]) * sy;
    position_array[v3 + 2] = center[2] + (position_array[v3 + 2] - center[2]) * sz;
  }
}

template <class T>
inline void GetRange(T* position_array, int num, T& min_x, T& min_y, T& min_z, T& max_x, T& max_y, T& max_z)
{
  min_x = min_y = min_z = 1e10;
  max_x = max_y = max_z = -1e10;
  for (int i = 0, v3 = 0; i < num; ++i, v3 += 3) {
    min_x = dj::Min(min_x, position_array[v3 + 0]);
    max_x = dj::Max(max_x, position_array[v3 + 0]);
    min_y = dj::Min(min_y, position_array[v3 + 1]);
    max_y = dj::Max(max_y, position_array[v3 + 1]);
    min_z = dj::Min(min_z, position_array[v3 + 2]);
    max_z = dj::Max(max_z, position_array[v3 + 2]);
  }
}

template <class T>
inline void GetRange(T* position_array, int num, T* min, T* max)
{
  GetRange(position_array, num,
           min[0], min[1], min[2],
           max[0], max[1], max[2]);
}

template <class T>
inline T GetMaxRange(T* position_array, int num)
{
  float min[3], max[3];
  GetRange<T>(position_array, num, min[0], min[1], min[2], max[0], max[1], max[2]);
  T max_range = dj::Max(max[0] - min[0], max[1] - min[1]);
  max_range = dj::Max(max_range, max[2] - min[2]);
  return max_range;
}

template <class T>
inline void RotateAroundX(T* position_array, int num, T angle_in_degree)
{
  T angle_in_radian = dj::Degree2Radian(angle_in_degree);
  for (int i = 0, v3 = 0; i < num; i++, v3 += 3) {
    T ty = position_array[v3 + 1] * cos(angle_in_radian) - position_array[v3 + 2] * sin(angle_in_radian);
    T tz = position_array[v3 + 1] * sin(angle_in_radian) + position_array[v3 + 2] * cos(angle_in_radian);
    position_array[v3 + 1] = ty;
    position_array[v3 + 2] = tz;
  }
}

template <class T>
inline void RotateAroundY(T* position_array, int num, T angle_in_degree)
{
  T angle_in_radian = dj::Degree2Radian(angle_in_degree);
  for (int i = 0, v3 = 0; i < num; i++, v3 += 3) {
    T tx = position_array[v3 + 0] *  cos(angle_in_radian) + position_array[v3 + 2] * sin(angle_in_radian);
    T tz = position_array[v3 + 0] * -sin(angle_in_radian) + position_array[v3 + 2] * cos(angle_in_radian);
    position_array[v3 + 0] = tx;
    position_array[v3 + 2] = tz;
  }
}

template <class T>
inline void RotateAroundZ(T* position_array, int num, T angle_in_degree)
{
  T angle_in_radian = dj::Degree2Radian(angle_in_degree);
  for (int i = 0, v3 = 0; i < num; i++, v3 += 3) {
    T tx = position_array[v3 + 0] * cos(angle_in_radian) - position_array[v3 + 1] * sin(angle_in_radian);
    T ty = position_array[v3 + 0] * sin(angle_in_radian) + position_array[v3 + 1] * cos(angle_in_radian);
    position_array[v3 + 0] = tx;
    position_array[v3 + 1] = ty;
  }
}

/// adapted from http://geomalgorithms.com/a07-_distance.html
template <class T>
inline T Segment2SegmentDistance(T* e00, T* e01, T* e10, T* e11, T& r, T&s, T* N = NULL)
{
  const T SMALL_NUM = T(1e-8);
  T u[3], v[3], w[3];
  // Vector   u = S1.P1 - S1.P0;
  u[0] = e01[0] - e00[0];
  u[1] = e01[1] - e00[1];
  u[2] = e01[2] - e00[2];
  // Vector   v = S2.P1 - S2.P0;
  v[0] = e11[0] - e10[0];
  v[1] = e11[1] - e10[1];
  v[2] = e11[2] - e10[2];

  // Vector   w = S1.P0 - S2.P0;
  w[0] = e00[0] - e10[0];
  w[1] = e00[1] - e10[1];
  w[2] = e00[2] - e10[2];
  T a = dj::Dot3(u, u);        // always >= 0
  T b = dj::Dot3(u, v);
  T c = dj::Dot3(v, v);        // always >= 0
  T d = dj::Dot3(u, w);
  T e = dj::Dot3(v, w);
  T D = a * c - b * b;    // always >= 0
  T sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
  T tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

  // compute the line parameters of the two closest points
  if (D < SMALL_NUM) { // the lines are almost parallel
    sN = T(0.0);         // force using point P0 on segment S1
    sD = T(1.0);         // to prevent possible division by 0.0 later
    tN = e;
    tD = c;
  } else {               // get the closest points on the infinite lines
    sN = (b * e - c * d);
    tN = (a * e - b * d);
    if (sN < T(0.0)) {        // sc < 0 => the s=0 edge is visible
      sN = T(0.0);
      tN = e;
      tD = c;
    } else if (sN > sD) { // sc > 1  => the s=1 edge is visible
      sN = sD;
      tN = e + b;
      tD = c;
    }
  }

  if (tN < T(0.0)) {            // tc < 0 => the t=0 edge is visible
    tN = T(0.0);
    // recompute sc for this edge
    if (-d < T(0.0))
      sN = T(0.0);
    else if (-d > a)
      sN = sD;
    else {
      sN = -d;
      sD = a;
    }
  } else if (tN > tD) {    // tc > 1  => the t=1 edge is visible
    tN = tD;
    // recompute sc for this edge
    if ((-d + b) < T(0.0))
      sN = T(0);
    else if ((-d + b) > a)
      sN = sD;
    else {
      sN = (-d +  b);
      sD = a;
    }
  }
  // finally do the division to get sc and tc
  sc = (dj::Abs(sN) < SMALL_NUM ? T(0.0) : sN / sD);
  tc = (dj::Abs(tN) < SMALL_NUM ? T(0.0) : tN / tD);

  // get the difference of the two closest points
  // Vector   dP = w + (sc * u) - (tc * v);  // =  S1(sc) - S2(tc)
  T dP[3];
  dP[0] = w[0] + sc * u[0] - tc * v[0];
  dP[1] = w[1] + sc * u[1] - tc * v[1];
  dP[2] = w[2] + sc * u[2] - tc * v[2];
  r = sc;
  s = tc;

  if (N != NULL) {
    N[0] = dP[0];
    N[1] = dP[1];
    N[2] = dP[2];
  }
  return dP[0] * dP[0] + dP[1] * dP[1] + dP[2] * dP[2];
}

inline
float PointToTriangleDistance2(const float *p, float *p0, float *p1, float *p2, float *pr, float &b0, float &b1, float &b2, float &weight, float bound)
{
  float e[3], e0[3], e1[3], e2[3], n[3], temp[3];
  (void) e1;
  for (int i = 0; i < 3; i++) {
    e [i] = p [i] - p0[i];
    e0[i] = p1[i] - p0[i];
    e1[i] = p2[i] - p1[i];
    e2[i] = p0[i] - p2[i];
  }

  dj::Cross3(e0, e2, n);
  float area = dj::Normalize3(n);
  weight = -dj::Dot3(e, n);

  if (weight > fabsf(bound) )	{return 99999;}

  //if(tag && weight<0.5)	return ;

  for (int i = 0; i < 3; i++)
    pr[i] = p[i] + weight * n[i];

  float er0[3], er1[3], er2[3];
  for (int i = 0; i < 3; i++) {
    er0[i] = pr[i] - p0[i];
    er1[i] = pr[i] - p1[i];
    er2[i] = pr[i] - p2[i];
  }

  dj::Cross3(er1, er2, temp);
  b0 = -dj::Dot3(temp, n);
  dj::Cross3(er2, er0, temp);
  b1 = -dj::Dot3(temp, n);
  //Cross(er0, er1, temp);
  b2 = area - b0 - b1;

  if (b0 >= 0 && b1 >= 0 && b2 >= 0) {
    float inv_sum = 1 / (b0 + b1 + b2);
    b0 *= inv_sum;
    b1 *= inv_sum;
    b2 *= inv_sum;
  } else {
    float *pp0, *pp1, e[3];
    float *bb0, *bb1;
    if (b0 < 0) {
      b0 = 0;
      if (p1 < p2)	{pp0 = p1; pp1 = p2; bb0 = &b1; bb1 = &b2; }
      else		{pp0 = p2; pp1 = p1; bb0 = &b2; bb1 = &b1; }
    } else if (b1 < 0) {
      b1 = 0;
      if (p0 < p2)	{pp0 = p2; pp1 = p0; bb0 = &b2; bb1 = &b0; }
      else		{pp0 = p0; pp1 = p2; bb0 = &b0; bb1 = &b2; }
    } else if (b2 < 0) {
      b2 = 0;
      if (p0 < p1)	{pp0 = p0; pp1 = p1; bb0 = &b0; bb1 = &b1; }
      else		{pp0 = p1; pp1 = p0; bb0 = &b1; bb1 = &b0; }
    }

    e[0] = pp1[0] - pp0[0];
    e[1] = pp1[1] - pp0[1];
    e[2] = pp1[2] - pp0[2];

    float e_length = dj::Normalize3(e);

    float ep[3];
    ep[0] = p[0] - pp0[0];
    ep[1] = p[1] - pp0[1];
    ep[2] = p[2] - pp0[2];


    *bb1 = dj::Dot3(ep, e) / e_length;
    *bb0 = 1 - *bb1;

    if (*bb1 < 0)	{*bb1 = 0; *bb0 = 1;}
    if (*bb0 < 0)	{*bb1 = 1; *bb0 = 0;}
  }

  pr[0] = p0[0] * b0 + p1[0] * b1 + p2[0] * b2;
  pr[1] = p0[1] * b0 + p1[1] * b1 + p2[1] * b2;
  pr[2] = p0[2] * b0 + p1[2] * b1 + p2[2] * b2;

  if (weight > 0) return  sqrtf((p[0] - pr[0]) * (p[0] - pr[0]) + (p[1] - pr[1]) * (p[1] - pr[1]) + (p[2] - pr[2]) * (p[2] - pr[2]));
  else		 return -sqrtf((p[0] - pr[0]) * (p[0] - pr[0]) + (p[1] - pr[1]) * (p[1] - pr[1]) + (p[2] - pr[2]) * (p[2] - pr[2]));
}

inline void PrintVec3(const char* name, const float* vec)
{
  printf("%s = [%f, %f, %f]\n", name, vec[0], vec[1], vec[2]);
}

inline
float PointToTriangleSignDistanceWithNormal(const float *p, float *p0, float *p1, float *p2,
                                            float* normal, float* triangle_vector, float triangle_area,
                                            float *pr, float* bary,
                                            float& dist_to_triangle_plane, float threshold)
{
  float e[3] = {
    p[0] - p0[0],
    p[1] - p0[1],
    p[2] - p0[2]
  };
  dist_to_triangle_plane = dj::Dot3(e, normal);
  //  float abs_dist = dj::Abs(dist_to_triangle_plane);
  if (dist_to_triangle_plane > dj::Abs(threshold))
    //    return abs_dist;
    return 99999;

  for (int i = 0; i < 3; i++) {
    pr[i] = p[i] - dist_to_triangle_plane * normal[i];
  }

  float er0[3] = {
    pr[0] - p0[0],
    pr[1] - p0[1],
    pr[2] - p0[2]
  };
  float vec[2][3];
  dj::Cross3(triangle_vector, er0, vec[0]);
  dj::Cross3(er0, triangle_vector + 3, vec[1]);
  bary[1] = dj::Dot3(normal, vec[1]);
  bary[2] = dj::Dot3(normal, vec[0]);
  bary[0] = triangle_area - bary[1] - bary[2];;

  if (bary[0] >= 0 && bary[1] >= 0 && bary[2] >= 0) {
    float inv_area = 1.0f / triangle_area;
    bary[1] *= inv_area;
    bary[2] *= inv_area;
    bary[0] = 1.0f - bary[1] - bary[2];;
    ASSERT(bary[0] >= -1e-6f && bary[0] <= 1.000001f, P(bary[0]));
    ASSERT(bary[1] >= -1e-6f && bary[1] <= 1.000001f, P(bary[1]));
    ASSERT(bary[2] >= -1e-6f && bary[2] <= 1.000001f, P(bary[2]));
    return dist_to_triangle_plane;
  } else {
    float& b0 = bary[0];
    float& b1 = bary[1];
    float& b2 = bary[2];
    float *pp0 = NULL, *pp1 = NULL, e[3];
    float *bb0 = NULL, *bb1 = NULL;

    if (b0 < 0) {
      b0 = 0;
      pp0 = p1; pp1 = p2; bb0 = &b1; bb1 = &b2;
    } else if (b1 < 0) {
      b1 = 0;
      pp0 = p2; pp1 = p0; bb0 = &b2; bb1 = &b0;
    } else if (b2 < 0) {
      b2 = 0;
      pp0 = p0; pp1 = p1; bb0 = &b0; bb1 = &b1;
    } else {
      ASSERT(false, P(b0, b1, b2));
    }

    e[0] = pp1[0] - pp0[0];
    e[1] = pp1[1] - pp0[1];
    e[2] = pp1[2] - pp0[2];

    float e_length_square = dj::Dot3(e, e);//dj::Normalize3(e);

    float ep[3];
    ep[0] = p[0] - pp0[0];
    ep[1] = p[1] - pp0[1];
    ep[2] = p[2] - pp0[2];


    *bb1 = dj::Dot3(ep, e) / e_length_square;
    *bb0 = 1 - *bb1;

    if (*bb1 < 0)	{*bb1 = 0; *bb0 = 1;}
    if (*bb0 < 0)	{*bb1 = 1; *bb0 = 0;}
    pr[0] = p0[0] * b0 + p1[0] * b1 + p2[0] * b2;
    pr[1] = p0[1] * b0 + p1[1] * b1 + p2[1] * b2;
    pr[2] = p0[2] * b0 + p1[2] * b1 + p2[2] * b2;

    float new_distance = sqrtf((p[0] - pr[0]) * (p[0] - pr[0]) + (p[1] - pr[1]) * (p[1] - pr[1]) + (p[2] - pr[2]) * (p[2] - pr[2]));
    return (dist_to_triangle_plane > 0) ? new_distance : -new_distance;
  }
}

inline
float PointToTriangleDistanceWithNormal(const float *p,
                                        float** tri_pos,
                                        float* normal, float* triangle_vector, float triangle_area,
                                        float *pr, float* bary, float threshold)
{
  float e[3] = {
    p[0] - tri_pos[0][0],
    p[1] - tri_pos[0][1],
    p[2] - tri_pos[0][2]
  };
  float dist_to_triangle_plane = dj::Dot3(e, normal);
  float abs_dist = dj::Abs(dist_to_triangle_plane);
  if (abs_dist > threshold) return abs_dist;

  for (int i = 0; i < 3; i++) {
    pr[i] = p[i] - dist_to_triangle_plane * normal[i];
  }

  float er0[3] = {
    pr[0] - tri_pos[0][0],
    pr[1] - tri_pos[0][1],
    pr[2] - tri_pos[0][2]
  };
  float vec[2][3];
  dj::Cross3(triangle_vector, er0, vec[0]);
  dj::Cross3(er0, triangle_vector + 3, vec[1]);
  bary[1] = dj::Dot3(normal, vec[1]);
  bary[2] = dj::Dot3(normal, vec[0]);
  bary[0] = triangle_area - bary[1] - bary[2];

  int i0, i1;
  if (bary[0] < 0) {
    bary[0] = 0;
    i0 = 1;
    i1 = 2;
  } else if (bary[1] < 0) {
    bary[1] = 0;
    i0 = 0;
    i1 = 2;
  } else if (bary[2] < 0) {
    bary[2] = 0;
    i0 = 0;
    i1 = 1;
  } else {
    ASSERT(bary[0] >= 0 && bary[1] >= 0 && bary[2] >= 0);
    float inv_area = 1.0f / triangle_area;
    bary[1] *= inv_area;
    bary[2] *= inv_area;
    bary[0] = 1.0f - bary[1] - bary[2];;
    ASSERT(bary[0] >= 0.0f && bary[0] <= 1.0f);
    ASSERT(bary[1] >= 0.0f && bary[1] <= 1.0f);
    ASSERT(bary[2] >= 0.0f && bary[2] <= 1.0f);
    return abs_dist;
  }

  e[0] = tri_pos[i1][0] - tri_pos[i0][0];
  e[1] = tri_pos[i1][1] - tri_pos[i0][1];
  e[2] = tri_pos[i1][2] - tri_pos[i0][2];
  float e_length_square = dj::Dot3(e, e);

  float ep[3] = {
    p[0] - tri_pos[i0][0],
    p[1] - tri_pos[i0][1],
    p[2] - tri_pos[i0][2]
  };

  bary[i1] = dj::Dot3(ep, e) / e_length_square;
  if (bary[i1] <= 0) {
    bary[i1] = 0;
    bary[i0] = 1;
    pr[0] = tri_pos[i0][0];
    pr[1] = tri_pos[i0][1];
    pr[2] = tri_pos[i0][2];
  } else if (bary[i1] >= 1) {
    bary[i1] = 1;
    bary[i0] = 0;
    pr[0] = tri_pos[i1][0];
    pr[1] = tri_pos[i1][1];
    pr[2] = tri_pos[i1][2];
  } else {
    bary[i0] = 1 - bary[i1];
    pr[0] = tri_pos[0][0] * bary[0] + tri_pos[1][0] * bary[1] + tri_pos[2][0] * bary[2];
    pr[1] = tri_pos[0][1] * bary[0] + tri_pos[1][1] * bary[1] + tri_pos[2][1] * bary[2];
    pr[2] = tri_pos[0][2] * bary[0] + tri_pos[1][2] * bary[1] + tri_pos[2][2] * bary[2];
  }
  float diff[3] = {
    p[0] - pr[0],
    p[1] - pr[1],
    p[2] - pr[2]
  };
  return sqrtf(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
}


inline
float PointToTriangleDistance(const float *p, float *p0, float *p1, float *p2, float *pr, float* bary)
{
  float& b0 = bary[0];
  float& b1 = bary[1];
  float& b2 = bary[2];
  float weight;
  float e[3], e0[3], e1[3], e2[3], n[3], temp[3];
  (void) e1;
  //  for (int i = 0; i < 3; i++) {
  //    ASSERT(p[i] == p[i]);
  //    ASSERT(p0[i] == p0[i]);
  //    ASSERT(p1[i] == p1[i]);
  //    ASSERT(p2[i] == p2[i]);
  //  }
  for (int i = 0; i < 3; i++) {
    e [i] = p [i] - p0[i];
    e0[i] = p1[i] - p0[i];
    e1[i] = p2[i] - p1[i];
    e2[i] = p0[i] - p2[i];
    //    ASSERT(e[i] == e[i]);
    //    ASSERT(e0[i] == e0[i]);
    //    ASSERT(e1[i] == e1[i]);
    //    ASSERT(e2[i] == e2[i]);
  }

  dj::Cross3(e0, e2, n);
  float area = dj::Normalize3(n);
  weight = -dj::Dot3(e, n);

  for (int i = 0; i < 3; i++)
    pr[i] = p[i] + weight * n[i];

  float er0[3], er1[3], er2[3];
  for (int i = 0; i < 3; i++) {
    er0[i] = pr[i] - p0[i];
    er1[i] = pr[i] - p1[i];
    er2[i] = pr[i] - p2[i];
  }

  dj::Cross3(er1, er2, temp);
  b0 = -dj::Dot3(temp, n);
  dj::Cross3(er2, er0, temp);
  b1 = -dj::Dot3(temp, n);
  b2 = area - b0 - b1;

  if (b0 >= 0 && b1 >= 0 && b2 >= 0) {
    float inv_sum = 1 / (b0 + b1 + b2);
    b0 *= inv_sum;
    b1 *= inv_sum;
    b2 *= inv_sum;
  } else {
    float *pp0 = NULL, *pp1 = NULL, e[3];
    float *bb0 = NULL, *bb1 = NULL;
    if (b0 < 0) {
      b0 = 0;
      if (p1 < p2)	{pp0 = p1; pp1 = p2; bb0 = &b1; bb1 = &b2; }
      else		{pp0 = p2; pp1 = p1; bb0 = &b2; bb1 = &b1; }
    } else if (b1 < 0) {
      b1 = 0;
      if (p0 < p2)	{pp0 = p2; pp1 = p0; bb0 = &b2; bb1 = &b0; }
      else		{pp0 = p0; pp1 = p2; bb0 = &b0; bb1 = &b2; }
    } else if (b2 < 0) {
      b2 = 0;
      if (p0 < p1)	{pp0 = p0; pp1 = p1; bb0 = &b0; bb1 = &b1; }
      else		{pp0 = p1; pp1 = p0; bb0 = &b1; bb1 = &b0; }
    }
    ASSERT(pp0 != NULL, P(b0, b1, b2));
    ASSERT(pp1 != NULL, P(b0, b1, b2));
    ASSERT(bb0 != NULL, P(b0, b1, b2));
    ASSERT(bb1 != NULL, P(b0, b1, b2));
    e[0] = pp1[0] - pp0[0];
    e[1] = pp1[1] - pp0[1];
    e[2] = pp1[2] - pp0[2];

    float e_length = dj::Normalize3(e);

    float ep[3];
    ep[0] = p[0] - pp0[0];
    ep[1] = p[1] - pp0[1];
    ep[2] = p[2] - pp0[2];


    *bb1 = dj::Dot3(ep, e) / e_length;
    *bb0 = 1 - *bb1;

    if (*bb1 < 0)	{*bb1 = 0; *bb0 = 1;}
    if (*bb0 < 0)	{*bb1 = 1; *bb0 = 0;}
  }

  pr[0] = p0[0] * b0 + p1[0] * b1 + p2[0] * b2;
  pr[1] = p0[1] * b0 + p1[1] * b1 + p2[1] * b2;
  pr[2] = p0[2] * b0 + p1[2] * b1 + p2[2] * b2;

  return sqrtf((p[0] - pr[0]) * (p[0] - pr[0]) + (p[1] - pr[1]) * (p[1] - pr[1]) + (p[2] - pr[2]) * (p[2] - pr[2]));
}

inline
float PointToTriangleDistance(const float *p, float *p0, float *p1, float *p2, float *pr, float &b0, float &b1, float &b2, float &weight, float bound)
{
  float e[3], e0[3], e1[3], e2[3], n[3], temp[3];
  (void) e1;
  for (int i = 0; i < 3; i++) {
    e [i] = p [i] - p0[i];
    e0[i] = p1[i] - p0[i];
    e1[i] = p2[i] - p1[i];
    e2[i] = p0[i] - p2[i];
  }

  dj::Cross3(e0, e2, n);
  //  printf("%f, %f, %f\n", n[0], n[1], n[2]);
  float area = dj::Normalize3(n);
  weight = -dj::Dot3(e, n);
  if (weight > fabsf(bound) )	{return 99999;}

  for (int i = 0; i < 3; i++)
    pr[i] = p[i] + weight * n[i];
  //  printf("%f, %f, %f\n", pr[0], pr[1], pr[2]);

  float er0[3], er1[3], er2[3];
  for (int i = 0; i < 3; i++) {
    er0[i] = pr[i] - p0[i];
    er1[i] = pr[i] - p1[i];
    er2[i] = pr[i] - p2[i];
  }

  dj::Cross3(er1, er2, temp);
  b0 = -dj::Dot3(temp, n);
  dj::Cross3(er2, er0, temp);
  b1 = -dj::Dot3(temp, n);
  //Cross(er0, er1, temp);
  b2 = area - b0 - b1;

  if (b0 >= 0 && b1 >= 0 && b2 >= 0) {
    float inv_sum = 1 / (b0 + b1 + b2);
    b0 *= inv_sum;
    b1 *= inv_sum;
    b2 *= inv_sum;
    return weight;
  } else {
    float *pp0 = NULL, *pp1 = NULL, e[3];
    float *bb0 = NULL, *bb1 = NULL;
    if (b0 < 0) {
      b0 = 0;
      if (p1 < p2)	{pp0 = p1; pp1 = p2; bb0 = &b1; bb1 = &b2; }
      else		{pp0 = p2; pp1 = p1; bb0 = &b2; bb1 = &b1; }
    } else if (b1 < 0) {
      b1 = 0;
      if (p0 < p2)	{pp0 = p2; pp1 = p0; bb0 = &b2; bb1 = &b0; }
      else		{pp0 = p0; pp1 = p2; bb0 = &b0; bb1 = &b2; }
    } else if (b2 < 0) {
      b2 = 0;
      if (p0 < p1)	{pp0 = p0; pp1 = p1; bb0 = &b0; bb1 = &b1; }
      else		{pp0 = p1; pp1 = p0; bb0 = &b1; bb1 = &b0; }
    }

    e[0] = pp1[0] - pp0[0];
    e[1] = pp1[1] - pp0[1];
    e[2] = pp1[2] - pp0[2];

    float e_length = dj::Normalize3(e);

    float ep[3];
    ep[0] = p[0] - pp0[0];
    ep[1] = p[1] - pp0[1];
    ep[2] = p[2] - pp0[2];


    *bb1 = dj::Dot3(ep, e) / e_length;
    *bb0 = 1 - *bb1;

    if (*bb1 < 0)	{*bb1 = 0; *bb0 = 1;}
    if (*bb0 < 0)	{*bb1 = 1; *bb0 = 0;}

    pr[0] = p0[0] * b0 + p1[0] * b1 + p2[0] * b2;
    pr[1] = p0[1] * b0 + p1[1] * b1 + p2[1] * b2;
    pr[2] = p0[2] * b0 + p1[2] * b1 + p2[2] * b2;
    if (weight > 0) return  sqrtf((p[0] - pr[0]) * (p[0] - pr[0]) + (p[1] - pr[1]) * (p[1] - pr[1]) + (p[2] - pr[2]) * (p[2] - pr[2]));
    else		 return -sqrtf((p[0] - pr[0]) * (p[0] - pr[0]) + (p[1] - pr[1]) * (p[1] - pr[1]) + (p[2] - pr[2]) * (p[2] - pr[2]));
  }
}

template <class T>
inline void EigenDecompose2(const T* m, T* eigen, T* eigen_vec1, T* eigen_vec2)
{
#if 1
  //  void solve_eig(double A, double B, double C, double D,
  //      double* lambda1, double *v1x, double*v1y,
  //      double* lambda2, double *v2x, double*v2y )
  //  {
  //      if(B*C <= tolerance  ) {
  //          *lambda1 = A; *v1x = 1; *v1y = 0;
  //          *lambda2 = D; *v2x = 0; *v2y = 1;
  //          return;
  //      }

  //      double tr = A + D;
  //      double det = A * D - B * C;
  //      double S = sqrt( square(tr/2) - det );
  //      *lambda1 = tr/2 + S;
  //      *lambda2 = tr/2 - S;

  //      double SS = sqrt(max(square((A-D)/2) + B * C, 0.0) );
  //      if( A - D < 0 ) {
  //          *v1x = C;
  //          *v1y = - (A-D)/2 + SS;
  //          *v2x = + (A-D)/2 - SS;
  //          *v2y = B;
  //      } else {
  //          *v2x = C;
  //          *v2y = - (A-D)/2 - SS;
  //          *v1x = + (A-D)/2 + SS;
  //          *v1y = B;
  //      }

  //      double n1 = sqrt(square(*v1x)+square(*v1y));
  //      *v1x /= n1; *v1y /= n1;
  //      double n2 = sqrt(square(*v2x)+square(*v2y));
  //      *v2x /= n2; *v2y /= n2;
  //  }

  const T kTolerance = T(1e-10);
  T reverse_diag_prod = m[1] * m[2];
  if (reverse_diag_prod <= kTolerance) {
    eigen[0] = m[0]; eigen_vec1[0] = 1; eigen_vec1[1] = 0;
    eigen[1] = m[3]; eigen_vec2[0] = 0; eigen_vec2[1] = 1;
    return;
  }

  T tr = m[0] + m[3]; //A + D;
  T det = m[0] * m[3] - reverse_diag_prod;
  T S = sqrt(dj::Square(tr / 2) - det);
  eigen[0] = tr / 2 + S;
  eigen[1] = tr / 2 - S;
  T SS = sqrt(dj::Max(dj::Square((m[0] - m[3]) / 2) + reverse_diag_prod, T(0.0)) );
  if (m[0] - m[3] < 0) {
    eigen_vec1[0] = m[2];
    eigen_vec1[1] = -(m[0] - m[3]) / 2 + SS;
    eigen_vec2[0] = +(m[0] - m[3]) / 2 - SS;
    eigen_vec2[1] = m[1];
  } else {
    eigen_vec2[0] = m[2];
    eigen_vec2[1] = -(m[0] - m[3]) / 2 - SS;
    eigen_vec1[0] = +(m[0] - m[3]) / 2 + SS;
    eigen_vec1[1] = m[1];
  }
  dj::Normalize2(eigen_vec1);
  dj::Normalize2(eigen_vec2);
#else
  T a = T(1.0);
  T b = -(m[3] + m[0]);
  T c = m[0] * m[3] - m[1] * m[2];
  T delta = b * b - 4 * a * c;
  if (delta < 0) {
    std::cerr << "EigenDecompose2() => FATAL: eigen decomposition failed with delta " << delta << "!!" << std::endl;
    exit(0);
  } else {
    delta = std::sqrt(delta);
    a *= T(2.0);
    a = T(1.0) / a;
    eigen[0] = (-b + delta) * a;
    eigen[1] = (-b - delta) * a;
    eigen_vec1[0] = eigen_vec2[0] = -b;
    eigen_vec1[1] = a - eigen[0];
    eigen_vec2[1] = a - eigen[1];
    dj::Normalize2(eigen_vec1);
    dj::Normalize2(eigen_vec2);
  }
#endif
}

template <class T, int n>
void gaussj(T (*a)[n], T (*b)[n])
{
  int i = 0, icol = 0, irow = 0, j = 0, k = 0, l = 0, ll = 0;
  T big, dum, pivinv;
  int indxc[n];
  int indxr[n];
  int ipiv[n];

  for (j = 0; j < n; j++) ipiv[j] = 0;
  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 1)
        for (k = 0; k < n; k++) {
          if (ipiv[k] == 0) {
            if (abs(a[j][k]) >= big) {
              big = abs(a[j][k]);
              irow = j;
              icol = k;
            }
          }
        }
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l = 0; l < n; l++) dj::Swap(a[irow][l], a[icol][l]);
      for (l = 0; l < n; l++) dj::Swap(b[irow][l], b[icol][l]);
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (a[icol][icol] == 0.0) printf("gaussj: Singular Matrix\n");
    pivinv = 1.0 / a[icol][icol];
    a[icol][icol] = 1.0;
    for (l = 0; l < n; l++) a[icol][l] *= pivinv;
    for (l = 0; l < n; l++) b[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
        dum = a[ll][icol];
        a[ll][icol] = 0.0;
        for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
        for (l = 0; l < n; l++) b[ll][l] -= b[icol][l] * dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
        dj::Swap(a[k][indxr[l]], a[k][indxc[l]]);
  }
}

template <class T, int n>
void GaussianElimination(T *a, T *b)
{
  int indxc[n];
  int indxr[n];
  int ipiv[n];
  int i, icol, irow, j, k, l, ll;
  T big, dum, pivinv, temp;

  for (j = 0; j < n; j++) ipiv[j] = 0;
  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 1)
        for (k = 0; k < n; k++) {
          if (ipiv[k] == 0) {
            if (fabs(a[j * n + k]) >= big) {
              big = fabs(a[j * n + k]);
              irow = j;
              icol = k;
            }
          }
        }
    ++(ipiv[icol]);

    if (irow != icol) {
      for (l = 0; l < n; l++) {temp = a[irow * n + l]; a[irow * n + l] = a[icol * n + l]; a[icol * n + l] = temp;}
      temp = b[irow]; b[irow] = b[icol]; b[icol] = temp;
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (a[icol * n + icol] == 0.0) printf("Error: Singular Matrix in Gaussian_Elimination.");
    pivinv = 1.0 / a[icol * n + icol];
    a[icol * n + icol] = 1.0;
    for (l = 0; l < n; l++) a[icol * n + l] *= pivinv;
    b[icol] *= pivinv;

    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
        dum = a[ll * n + icol];
        a[ll * n + icol] = 0.0;
        for (l = 0; l < n; l++) a[ll * n + l] -= a[icol * n + l] * dum;
        b[ll] -= b[icol] * dum;
      }
  }

  for (l = n - 1; l > 1; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++) {
        temp = a[k * n + indxr[l]];
        a[k * n + indxr[l]] = a[k * n + indxc[l]];
        a[k * n + indxc[l]] = temp;
      }
  }
}

template <class T, int n>
void GaussianElimination(T (*a)[n], T *b) {
  GaussianElimination<T, n>((T*) a, b);
}

template <class T>
inline void Rotate(T angele_in_radian, const T* rotation_axis, T* matrix)
{
  T (*m)[3] = (T (*)[3]) matrix;
  T sin = std::sin(angele_in_radian);
  T cos = std::cos(angele_in_radian);
  m[0][0] = cos + rotation_axis[0] * rotation_axis[0] * (1 - cos);
  m[0][1] = rotation_axis[0] * rotation_axis[1] * (1 - cos) - rotation_axis[2] * sin;
  m[0][2] = rotation_axis[0] * rotation_axis[2] * (1 - cos) + rotation_axis[1] * sin;
  m[1][0] = rotation_axis[0] * rotation_axis[1] * (1 - cos) + rotation_axis[2] * sin;
  m[1][1] = cos + rotation_axis[1] * rotation_axis[1] * (1 - cos);
  m[1][2] = rotation_axis[1] * rotation_axis[2] * (1 - cos) - rotation_axis[0] * sin;
  m[2][0] = rotation_axis[0] * rotation_axis[2] * (1 - cos) - rotation_axis[1] * sin;
  m[2][1] = rotation_axis[1] * rotation_axis[2] * (1 - cos) + rotation_axis[0] * sin;
  m[2][2] = cos + rotation_axis[2] * rotation_axis[2] * (1 - cos);
}

//template <class T, int d>
//std::ostream& operator<<(std::ostream& out, dj::Vec<T, d> vec);
typedef std::vector<int> vectori;
typedef std::vector<unsigned int> vectorui;
typedef std::vector<float> vectorf;
typedef std::vector<double> vectord;
typedef std::vector<bool> vectorb;
typedef std::vector<char> vectorc;
typedef std::vector<unsigned char> vectoruc;


typedef FixedMatrixWrapper<unsigned int, 2, 2> Wmat2ui;
typedef FixedMatrixWrapper<int, 2, 2> Wmat2i;
typedef FixedMatrixWrapper<float, 2, 2> Wmat2f;
typedef FixedMatrixWrapper<double, 2, 2> Wmat2d;

typedef FixedMatrixWrapper<unsigned int, 3, 3> Wmat3ui;
typedef FixedMatrixWrapper<int, 3, 3> Wmat3i;
typedef FixedMatrixWrapper<float, 3, 3> Wmat3f;
typedef FixedMatrixWrapper<double, 3, 3> Wmat3d;

typedef FixedMatrixWrapper<unsigned int, 4, 4> Wmat4ui;
typedef FixedMatrixWrapper<int, 4, 4> Wmat4i;
typedef FixedMatrixWrapper<float, 4, 4> Wmat4f;
typedef FixedMatrixWrapper<double, 4, 4> Wmat4d;

typedef FixedMatrix<unsigned int, 2, 2> Mat2ui;
typedef FixedMatrix<int, 2, 2> Mat2i;
typedef FixedMatrix<float, 2, 2> Mat2f;
typedef FixedMatrix<double, 2, 2> Mat2d;

typedef FixedMatrix<unsigned int, 3, 3> Mat3ui;
typedef FixedMatrix<int, 3, 3> Mat3i;
typedef FixedMatrix<float, 3, 3> Mat3f;
typedef FixedMatrix<double, 3, 3> Mat3d;

typedef FixedMatrix<unsigned int, 4, 4> Mat4ui;
typedef FixedMatrix<int, 4, 4> Mat4i;
typedef FixedMatrix<float, 4, 4> Mat4f;
typedef FixedMatrix<double, 4, 4> Mat4d;
}
#endif // VECOTR_LIB_H_
