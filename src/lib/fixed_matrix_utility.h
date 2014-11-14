#ifndef FIXED_MATRIX_UTILITY_H_
#define FIXED_MATRIX_UTILITY_H_
#include "macro_constant.h"
#include "utility_function.h"
#include "fixed_vector_utility.h"
namespace dj
{

template <class T, int row_num, int col_num, class T1, class T2>
inline void Transpose(const T1* matrix_in, T2* transpose_out)
{
  if ((void*) matrix_in != (void*) transpose_out) {
    T (*matrix)[col_num] = (T (*)[col_num]) matrix_in;
    T (*transpose)[row_num] = (T (*)[row_num]) transpose_out;
    for (int row = 0; row < row_num; ++row) {
      for (int col = 0; col < col_num; ++col) {
        transpose[row][col] = matrix[col][row];
      }
    }
  }
}

template <class T, int row_num, int col_num, class T1>
inline void Transpose(const T1* matrix)
{
  static_assert(row_num == col_num, "row num is not equal to col_num in Transpose()");
  T (*real_matrix)[col_num] = (T (*)[col_num]) matrix;
  for (int row = 0; row < row_num; ++row) {
    for (int col = row + 1; col < col_num; ++col) {
      T tmp = real_matrix[row][col];
      real_matrix[row][col] = real_matrix[col][row];
      real_matrix[col][row] = tmp;
    }
  }
}

template <class T, class T1>
inline void Transpose3x3(T1* matrix)
{
  T (*real_matrix)[3] = (T (*)[3]) matrix;
  Swap(real_matrix[0][1], real_matrix[1][0]);
  Swap(real_matrix[0][2], real_matrix[2][0]);
  Swap(real_matrix[1][2], real_matrix[2][1]);
}

template <class T>
inline T Determinant2(T* matrix)
{
  return matrix[0] * matrix[3] - matrix[1] * matrix[2];
}

template <class T>
inline void InverseMatrix2(T* matrix, T* inverse_matrix)
{
  T inv_det = T(1.0) / Determinant2<T>(matrix);
  inverse_matrix[0] = matrix[3] * inv_det;
  inverse_matrix[1] = matrix[1] * -inv_det;
  inverse_matrix[2] = matrix[2] * -inv_det;
  inverse_matrix[3] = matrix[0] * inv_det;
}

template <class T>
inline void MulVecMatrix3x3(const T* vec, const T (*matrix)[3], T* result)
{
  result[0] = matrix[0][0] * vec[0] + matrix[1][0] * vec[1] + matrix[2][0] * vec[2];
  result[1] = matrix[0][1] * vec[0] + matrix[1][1] * vec[1] + matrix[2][1] * vec[2];
  result[2] = matrix[0][2] * vec[0] + matrix[1][2] * vec[1] + matrix[2][2] * vec[2];
}

template <class T>
inline void MulMatrix3x3Vec(const T (*matrix)[3], const T* vec, T* result)
{
  result[0] =  matrix[0][0] * vec[0] + matrix[0][1] * vec[1] + matrix[0][2] * vec[2];
  result[1] =  matrix[1][0] * vec[0] + matrix[1][1] * vec[1] + matrix[1][2] * vec[2];
  result[2] =  matrix[2][0] * vec[0] + matrix[2][1] * vec[1] + matrix[2][2] * vec[2];
}

// result = left^T * right
template <class T, class T1, class T2, class T3>
inline void MulMatrixLeftTransposed3x3(const T1* left, const T2* right, T3* result)
{
  T (*a)[3] = (T (*)[3]) left;
  T (*b)[3] = (T (*)[3]) right;
  T (*real_result)[3] = (T (*)[3]) result;
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      real_result[row][col] = a[0][row] * b[0][col] + a[1][row] * b[1][col] + a[2][row] * b[2][col];
    }
  }
}

// result = left * right^T
template <class T, class T1, class T2, class T3>
inline void MulMatrixRightTransposed3x3(const T1* left, const T2* right, T3* result)
{
  T (*a)[3] = (T (*)[3]) left;
  T (*b)[3] = (T (*)[3]) right;
  T (*real_result)[3] = (T (*)[3]) result;
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      real_result[row][col] =  a[row][0] * b[col][0] + a[row][1] * b[col][1] + a[row][2] * b[col][2];
    }
  }
}

template <class T, class T1, class T2, class T3>
inline void MulMatrix3x3(const T1* left, const T2* right, T3* result)
{
  T (*a)[3] = (T (*)[3]) left;
  T (*b)[3] = (T (*)[3]) right;
  T (*real_result)[3] = (T (*)[3]) result;
  int row = 0, col = 0;
  real_result[row][col] =  a[row][0] * b[0][col] + a[row][1] * b[1][col] + a[row][2] * b[2][col];
  row = 0, col = 1;
  real_result[row][col] =  a[row][0] * b[0][col] + a[row][1] * b[1][col] + a[row][2] * b[2][col];
  row = 0, col = 2;
  real_result[row][col] =  a[row][0] * b[0][col] + a[row][1] * b[1][col] + a[row][2] * b[2][col];
  row = 1, col = 0;
  real_result[row][col] =  a[row][0] * b[0][col] + a[row][1] * b[1][col] + a[row][2] * b[2][col];
  row = 1, col = 1;
  real_result[row][col] =  a[row][0] * b[0][col] + a[row][1] * b[1][col] + a[row][2] * b[2][col];
  row = 1, col = 2;
  real_result[row][col] =  a[row][0] * b[0][col] + a[row][1] * b[1][col] + a[row][2] * b[2][col];
  row = 2, col = 0;
  real_result[row][col] =  a[row][0] * b[0][col] + a[row][1] * b[1][col] + a[row][2] * b[2][col];
  row = 2, col = 1;
  real_result[row][col] =  a[row][0] * b[0][col] + a[row][1] * b[1][col] + a[row][2] * b[2][col];
  row = 2, col = 2;
  real_result[row][col] =  a[row][0] * b[0][col] + a[row][1] * b[1][col] + a[row][2] * b[2][col];
}


template <class T, int row_num1, int col_num1, int col_num2, class T1, class T2, class T3>
inline void MulMatrix(T1* m1, T2* m2, T3* result)
{
  T (*matrix1)[col_num1] = (T (*)[col_num1]) m1;
  T (*matrix2)[col_num2] = (T (*)[col_num2]) m2;
  T (*result_matrix)[col_num2] = (T (*)[col_num2]) result;
  for (int row = 0; row < row_num1; ++row) {
    for (int col = 0; col < col_num2; ++col) {
      result_matrix[row][col] = 0;
      for (int i = 0; i < col_num1; ++i) {
        result_matrix[row][col] += matrix1[row][i] * matrix2[i][col];
      }
    }
  }
}

template <class T, int row_num, class T1, class T2, class T3>
inline void MulSquareMatrix(T1* m1, T2* m2, T3* result)
{
  MulMatrix<T, row_num, row_num, row_num>(m1, m2, result);
}

template <class T, int row_num, int col_num, class T1>
inline void MulMatrixVec(T1* m, const T* vec, T* result)
{
  T (*matrix)[col_num] = (T (*)[col_num]) m;
  for (int row = 0; row < row_num; ++row) {
    result[row] = 0;
    for (int col = 0; col < col_num; ++col) {
      result[row] += matrix[row][col] * vec[col];
    }
  }
}

template <class T, int row_num, int col_num, class T1, class T2, class T3>
inline void AddMatrix(T1* m1, T2* m2, T3* result)
{
  const int size = row_num * col_num;
  T* array1 = (T*) m1;
  T* array2 = (T*) m2;
  T* result_array = (T*) result;
  for (int i = 0; i < size; ++i) {
    result_array[i] = m1[i] + m2[i];
  }
}

template <class T, int row_num, int col_num, class T1, class T2, class T3>
inline void SubMatrix(T1* m1, T2* m2, T3* result)
{
  const int size = row_num * col_num;
  T* array1 = (T*) m1;
  T* array2 = (T*) m2;
  T* result_array = (T*) result;
  for (int i = 0; i < size; ++i) {
    result_array[i] = m1[i] - m2[i];
  }
}

template <class T>
inline T OneNorm(const T* A)
{
  T norm = T(0);
  for (int i = 0; i < 3; i++) {
    T columnAbsSum = Abs(A[i + 0]) + Abs(A[i + 3]) + Abs(A[i + 6]);
    if (columnAbsSum > norm)
      norm = columnAbsSum;
  }
  return norm;
}

template <class T>
inline T InfNorm(const T* A)
{
  T norm = 0.0;
  for (int i = 0; i < 3; i++) {
    T rowSum = Abs(A[3 * i + 0]) + Abs(A[3 * i + 1]) + Abs(A[3 * i + 2]);
    if (rowSum > norm)
      norm = rowSum;
  }
  return norm;
}

template <class T>
inline T PolarDecompose(const T* M, T * Q, T * S, T tolerance)
{
  T Mk[9];
  T Ek[9];
  T det, M_oneNorm, M_infNorm, E_oneNorm;

  // Mk = M^T
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Mk[3 * i + j] = M[3 * j + i];

  M_oneNorm = OneNorm(Mk);
  M_infNorm = InfNorm(Mk);

  do {
    T MadjTk[9];
    // row 2 x row 3
    Cross3(&(Mk[3]), &(Mk[6]), &(MadjTk[0]));
    // row 3 x row 1
    Cross3(&(Mk[6]), &(Mk[0]), &(MadjTk[3]));
    // row 1 x row 2
    Cross3(&(Mk[0]), &(Mk[3]), &(MadjTk[6]));

    det = Mk[0] * MadjTk[0] + Mk[1] * MadjTk[1] + Mk[2] * MadjTk[2];
    if (det > T(-1.0e-6) && det < T(1.0e-6)) {
      return det;
    }

    T MadjT_one = OneNorm(MadjTk);
    T MadjT_inf = InfNorm(MadjTk);

    T gamma = sqrt(sqrt((MadjT_one * MadjT_inf) / (M_oneNorm * M_infNorm)) / Abs(det));
    T g1 = gamma * T(0.5);
    T g2 = T(0.5) / (gamma * det);

    for (int i = 0; i < 9; i++) {
      Ek[i] = Mk[i];
      Mk[i] = g1 * Mk[i] + g2 * MadjTk[i];
      Ek[i] -= Mk[i];
    }

    E_oneNorm = OneNorm(Ek);
    M_oneNorm = OneNorm(Mk);
    M_infNorm = InfNorm(Mk);
  } while ( E_oneNorm > M_oneNorm * tolerance );

  // Q = Mk^T
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Q[3 * i + j] = Mk[3 * j + i];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      S[3 * i + j] = 0;
      for (int k = 0; k < 3; k++)
        S[3 * i + j] += Mk[3 * i + k] * M[3 * k + j];
    }

  // S must be symmetric; enforce the symmetry
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
      S[3 * i + j] = S[3 * j + i] = T(0.5) * (S[3 * i + j] + S[3 * j + i]);

  return (det);
}


template <class T, int row, int col>
struct FixedMatOp {
  static const int kSize = row * col;

  template <class T1>
  inline static void MulElement(const T* m1, const T1& mul, typename MulType<T, T1>::Type* result)
  {
    typedef typename MulType<T, T1>::Type ResultType;
    const T (*x)[col] = (T (*)[col]) m1;
    ResultType (*z)[col] = (ResultType (*)[col]) result;
    for (int i = 0; i < row; ++i) {
      for (int j = 0; j < col; ++j) {
        z[i][j] = mul * x[i][j];
      }
    }
  }

  template <class T1, int col2>
  inline static void Mul(const T* m1, const T1* m2, typename MulType<T, T1>::Type* result)
  {
    typedef typename MulType<T, T1>::Type ResultType;
    const T (*x)[col] = (T (*)[col]) m1;
    const T1 (*y)[col2] = (T1 (*)[col2]) m2;
    ResultType (*z)[col2] = (ResultType (*)[col2]) result;
    for (int i = 0; i < row; ++i) {
      for (int j = 0; j < col2; ++j) {
        z[i][j] = ResultType();
        for (int k = 0; k < col; ++k) {
          z[i][j] += x[i][k] * y[k][j];
        }
      }
    }
  }

  template <class T1>
  inline static void MulLeft(const T1* vec, const T* m, typename MulType<T, T1>::Type* result)
  {
    typedef typename MulType<T, T1>::Type ResultType;
    const T (*x)[col] = (T (*)[col]) m;
    for (int k = 0; k < col; ++k) {
      result[k] = ResultType();
      for (int i = 0; i < row; ++i) {
        result[k] += vec[i] * x[i][k];
      }
    }
  }

  template <class T1>
  inline static void Mul(const T* m1, const T1* vec, typename MulType<T, T1>::Type* result)
  {
    typedef typename MulType<T, T1>::Type ResultType;
    const T (*x)[col] = (T (*)[col]) m1;
    for (int i = 0; i < row; ++i) {
      result[i] = ResultType();
      for (int k = 0; k < col; ++k) {
        result[i] += x[i][k] * vec[k];
      }
    }
  }

  inline static void Plus(const T* m1, const T* m2, const T* result)
  {
    for (int i = 0; i < kSize; ++i) {
      result[i] = m1[i] + m2[i];
    }
  }

  inline static void Minus(const T* m1, const T* m2, const T* result)
  {
    for (int i = 0; i < kSize; ++i) {
      result[i] = m1[i] - m2[i];
    }
  }

  inline static bool Equal(const T* x1, const T* x2)
  {
    for (int i = 0; i < kSize; ++i) {
      if (x1[i] != x2[i]) return false;
    }
    return true;
  }

  inline static void Copy(const T* from, T* to)
  {
    for (int i = 0; i < kSize; ++i) {
      to[i] = from[i];
    }
  }

  inline static void Scale(T* src, T scale, T* dest)
  {
    for (int i = 0; i < kSize; ++i) {
      dest[i] = scale * src[i];
    }
  }

  inline static void Minus(const T* x, const T* y, T* result)
  {
    for (int i = 0; i < kSize; ++i) {
      result[i] = x[i] - y[i];
    }
  }

  inline static void Plus(const T* x, const T* y, T* result)
  {
    for (int i = 0; i < kSize; ++i) {
      result[i] = x[i] + y[i];
    }
  }

  inline static void Fill(T* x, const T value)
  {
    for (int i = 0; i < kSize; ++i) {
      x[i] = value;
    }
  }

  inline static T MinAbs(const T* x)
  {
    T min = std::numeric_limits<T>::max();
    for (int i = 0; i < kSize; i++) {
      T abs = dj::Abs(x[i]);
      if (min > abs) {
        min = abs;
      }
    }
    return min;
  }

  inline static T MaxAbs(const T* x)
  {
    T max = T(EPSILON);
    for (int i = 0; i < kSize; i++) {
      T abs = dj::Abs(x[i]);
      if (max < abs) {
        max = abs;
      }
    }
    return max;
  }

  inline static T Min(const T* x)
  {
    T min = std::numeric_limits<T>::max();
    for (int i = 0; i < kSize; i++) {
      if (min > x[i]) {
        min = x[i];
      }
    }
    return min;
  }

  inline static T Max(const T* x)
  {
    T max = std::numeric_limits<T>::min();
    for (int i = 0; i < kSize; i++) {
      if (max < x[i]) {
        max = x[i];
      }
    }
    return max;
  }
};


} // namespace
#endif // FIXED_MATRIX_UTILITY_H_
