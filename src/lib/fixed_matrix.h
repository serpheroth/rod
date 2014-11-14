#ifndef FIXED_MATRIX_H_
#define FIXED_MATRIX_H_
#include <type_traits>
#include "fixed_matrix_utility.h"
#include "fixed_vector.h"
#include "print_macro.h"
namespace  dj
{

template <class T, int col_num>
inline Vec<T, col_num, false> GetRow(const T (*matrix)[col_num], int row)
{
  typedef Vec<T, col_num, false> RowVector;
  return RowVector(matrix[row]);
}

template <class T, int row_num, int col_num>
inline Vec<T, row_num, false> GetColumn(const T (*matrix)[col_num], int col)
{
  typedef Vec<T, row_num, false> ColumnVector;
  ColumnVector column;
  for (int i = 0; i < row_num; ++i) {
    column[i] = matrix[i][col];
  }
  return column;
}

template <class T, int row_num, int col_num>
class FixedMatrix;

template <class T, int row_num, int col_num>
class FixedMatrixWrapper
{
  static_assert(row_num >= 1 && col_num >= 1, "dimension of matrix must be at least 1");
private:
  typedef FixedMatrix<T, row_num, col_num> Mat;
  typedef FixedMatrixWrapper<T, row_num, col_num> Wrapper;
  typedef Vec<T, col_num, false> RowVector;
  typedef Vec<T, row_num, false> ColumnVector;
public:
  enum { kRowNum = row_num };
  enum { kColNum = col_num };
  enum { kSize = kRowNum * kColNum };

  template <class T1>
  FixedMatrixWrapper(T1* matrix) {
    matrix_ = (T (*)[col_num]) matrix;
  }

  FixedMatrixWrapper(const Wrapper& o) {
    this->matrix_ = o.matrix_;
  }

  const Wrapper& operator=(const FixedMatrixWrapper<T, row_num, col_num>& o) {
    this->matrix_ = o->matrix_;
    return *this;
  }

  T* operator()() const {
    return (T*) matrix_;
  }

  RowVector operator[](int row) const {
    return GetRow<T, col_num>(matrix_, row);
  }
  RowVector Row(int row) const {
    return GetRow<T, col_num>(matrix_, row);
  }

  ColumnVector operator()(int col) const {
    return GetColumn<T, row_num, col_num>(matrix_, col);
  }

  ColumnVector Col(int col) const {
    return GetColumn<T, row_num, col_num>(matrix_, col);
  }

  T& operator()(int x, int y) const {
    return matrix_[x][y];
  }

  Mat operator+(const Wrapper& other) const {
    Mat result;
    FixedMatOp<T, row_num, col_num>::Plus((T*) matrix_, (T*) other.matrix_, result());
    return result;
  }

  Wrapper& operator+=(const Wrapper& other) {
    FixedMatOp<T, row_num, col_num>::Plus((T*) matrix_, (T*) other.matrix_, (T*) matrix_);
    return *this;
  }

  Wrapper& operator+=(const Mat& other) {
    FixedMatOp<T, row_num, col_num>::Plus((T*) matrix_, (T*) other.matrix_, (T*) matrix_);
    return *this;
  }

  Mat operator+(const Mat& other) const {
    Mat result;
    FixedMatOp<T, row_num, col_num>::Plus((T*) matrix_, (T*) other.matrix_, result());
    return result;
  }

  Mat operator-(const Mat& other) const {
    Mat result;
    FixedMatOp<T, row_num, col_num>::Minus((T*) matrix_, (T*) other.matrix_, result());
    return result;
  }
  Mat operator-(const Wrapper& other) const {
    Mat result;
    FixedMatOp<T, row_num, col_num>::Minus((T*) matrix_, (T*) other.matrix_, result());
    return result;
  }

  Wrapper& operator-=(const Wrapper& other) {
    FixedMatOp<T, row_num, col_num>::Minus((T*) matrix_, (T*) other.matrix_, (T*) matrix_);
    return *this;
  }

  Wrapper& operator-=(const Mat& other) {
    FixedMatOp<T, row_num, col_num>::Minus((T*) matrix_, (T*) other.matrix_, (T*) matrix_);
    return *this;
  }

  Mat operator/(const T& div) const {
    Mat result;
    FixedMatOp<T, row_num, col_num>::template MulElement<T>((T*) matrix_, T(1) / div, result());
    return result;
  }

  Wrapper& operator/=(const T& div) const {
    FixedMatOp<T, row_num, col_num>::template MulElement<T>((T*) matrix_, T(1) / div, (T*) matrix_);
    return *this;
  }


  Wrapper& Transpose(Wrapper& o) {
    dj::Transpose<T, row_num, col_num>(matrix_, o.matrix_);
    return o;
  }

  Wrapper& Transpose() {
    dj::Transpose<T, row_num, col_num>(matrix_);
    return *this;
  }

  template <class T1, int new_col_num>
  FixedMatrix<typename MulType<T, T1>::Type, row_num, new_col_num> operator*(const FixedMatrix<T1, col_num, new_col_num>& other) const {
    typedef typename MulType<T, T1>::Type ElementType;
    typedef FixedMatrix<ElementType, row_num, new_col_num> ResultMatrix;
    ResultMatrix result;
    FixedMatOp<T, row_num, col_num>::template Mul<T1, new_col_num>((T*) matrix_, (T1*) other.matrix_, (ElementType*) result.matrix_);
    return result;
  }

  template <class T1, int new_col_num>
  FixedMatrix<typename MulType<T, T1>::Type, row_num, new_col_num> operator*(const FixedMatrixWrapper<T1, col_num, new_col_num>& other) const {
    typedef typename MulType<T, T1>::Type ElementType;
    typedef FixedMatrix<ElementType, row_num, new_col_num> ResultMatrix;
    ResultMatrix result;
    FixedMatOp<T, row_num, col_num>::template Mul<T1, new_col_num>((T*) matrix_, other(), (ElementType*) result.matrix_);
    return result;
  }

  template <class T1>
  Vec<typename MulType<T, T1>::Type, row_num, false> operator*(T1* vec) const {
    typedef typename MulType<T, T1>::Type ElementType;
    typedef Vec<ElementType, row_num, false> ResultVector;
    ResultVector result;
    FixedMatOp<T, row_num, col_num>::template Mul<T1>((T*) matrix_, vec, result());
    return result;
  }

  template <class T1, bool wrapped>
  Vec<typename MulType<T, T1>::Type, row_num, false> operator*(const Vec<T1, col_num, wrapped>& vec) const {
    typedef typename MulType<T, T1>::Type ElementType;
    typedef Vec<ElementType, row_num, false> ResultVector;
    ResultVector result;
    FixedMatOp<T, row_num, col_num>::template Mul<T1>((T*) matrix_, vec(), result());
    return result;
  }

  template <class T1>
  Wrapper& operator*=(const T1& mul) {
    static_assert(std::is_same<T, typename MulType<T, T1>::Type>::value,
                  "Result type of multiplication must be the same as element type");
    FixedMatOp<T, row_num, col_num>::template MulElement<T1>((T*) matrix_, mul, (T*) matrix_);
    return *this;
  }

  T (*matrix_)[col_num];
private:
  FixedMatrixWrapper() {}
};


enum MatrixType {
  kZero,
  kIdentity
};

template <class T, int row_num, int col_num>
class FixedMatrix
{
  static_assert(row_num >= 1 && col_num >= 1, "dimension of matrix must be at least 1");
public:
  static const int kRowNum = row_num;
  static const int kColNum = col_num;
  static const int kSize = row_num * col_num;
private:
  typedef T (*Array)[col_num];
  typedef FixedMatrix<T, row_num, col_num> Mat;
  typedef FixedMatrixWrapper<T, row_num, col_num> Wrapper;
  typedef Vec<T, col_num, false> RowVector;
  typedef Vec<T, row_num, false> ColumnVector;
public:
  FixedMatrix(MatrixType type = kZero) {
    memset(matrix_, 0, sizeof(T) * row_num * col_num);
    if (type == kIdentity && row_num == col_num) {
      for (int i = 0; i < row_num; ++i) {
        matrix_[i][i] = T(1);
      }
    }
  }
  /*
  FixedMatrix(const std::initializer_list<T>& list) {
    T* m = (T*) matrix_;
    int i = 0;
    for (auto iter = list.begin(); i < kSize; ++iter, ++i) {
      m[i] = *iter;
    }
  }
  */
  template <class T1>
  FixedMatrix(T1* matrix) {
    FixedMatOp<T, row_num, col_num>::Copy((T*) matrix, (T*) matrix_);
  }

  template <class T1>
  FixedMatrix(const T1* matrix) {
    FixedMatOp<T, row_num, col_num>::Copy((T*) matrix, (T*) matrix_);
  }

  FixedMatrix(const Mat& o) {
    FixedMatOp<T, row_num, col_num>::Copy((T*) o.matrix_, (T*) matrix_);
  }

  FixedMatrix(const Wrapper& o) {
    FixedMatOp<T, row_num, col_num>::Copy(o.matrix_, matrix_);
  }

  const Mat& operator=(const Mat& o) {
    memcpy((void*) matrix_, (void*) o.matrix_, sizeof(T) * col_num * row_num);
    return *this;
  }

  T* operator()() const {
    return (T*) matrix_;
  }

  RowVector operator[](int row) const {
    return GetRow<T, col_num>((Array) matrix_, row);
  }
  RowVector Row(int row) const {
    return GetRow<T, col_num>((Array) matrix_, row);
  }

  ColumnVector operator()(int col) const {
    return GetColumn<T, row_num, col_num>((Array) matrix_, col);
  }
  ColumnVector Col(int col) const {
    return GetColumn<T, row_num, col_num>((Array) matrix_, col);
  }

  const T& operator()(int x, int y) const {
    return matrix_[x][y];
  }

  T& operator()(int x, int y) {
    return matrix_[x][y];
  }

  template <class MatrixTemplate>
  Mat operator+(const MatrixTemplate& o) const {
    Mat result;
    AddMatrix<T, row_num, col_num>(matrix_, o(), result.matrix_);
    return result;
  }

  template <class MatrixTemplate>
  Mat& operator+=(const MatrixTemplate& o) {
    AddMatrix<T, row_num, col_num>(matrix_, o(), matrix_);
    return *this;
  }

  template <class MatrixTemplate>
  Mat& Transpose(MatrixTemplate& o) const {
    dj::Transpose<T, row_num, col_num>(matrix_, o());
    return o;
  }

  Mat& Transpose() {
    dj::Transpose<T, row_num, col_num>(matrix_);
    return *this;
  }

  Mat operator+(const Wrapper& other) const {
    Mat result;
    FixedMatOp<T, row_num, col_num>::Plus((T*) matrix_, (T*) other.matrix_, result());
    return result;
  }

  Mat& operator+=(const Wrapper& other) {
    FixedMatOp<T, row_num, col_num>::Plus((T*) matrix_, (T*) other.matrix_, (T*) matrix_);
    return *this;
  }

  Mat& operator+=(const Mat& other) {
    FixedMatOp<T, row_num, col_num>::Plus((T*) matrix_, (T*) other.matrix_, (T*) matrix_);
    return *this;
  }

  Mat operator+(const Mat& other) const {
    Mat result;
    FixedMatOp<T, row_num, col_num>::Plus((T*) matrix_, (T*) other.matrix_, result());
    return result;
  }

  Mat operator-(const Mat& other) const {
    Mat result;
    FixedMatOp<T, row_num, col_num>::Minus((T*) matrix_, (T*) other.matrix_, result());
    return result;
  }
  Mat operator-(const Wrapper& other) const {
    Mat result;
    FixedMatOp<T, row_num, col_num>::Minus((T*) matrix_, (T*) other.matrix_, result());
    return result;
  }

  Mat& operator-=(const Wrapper& other) {
    FixedMatOp<T, row_num, col_num>::Minus((T*) matrix_, (T*) other.matrix_, (T*) matrix_);
    return *this;
  }

  Mat& operator-=(const Mat& other) {
    FixedMatOp<T, row_num, col_num>::Minus((T*) matrix_, (T*) other.matrix_, (T*) matrix_);
    return *this;
  }

  Mat operator/(const T& div) const {
    Mat result;
    FixedMatOp<T, row_num, col_num>::template MulElement<T>((T*) matrix_, T(1) / div, result());
    return result;
  }

  Mat& operator/=(const T& div) const {
    FixedMatOp<T, row_num, col_num>::template MulElement<T>((T*) matrix_, T(1) / div, (T*) matrix_);
    return *this;
  }

  template <class T1, int new_col_num>
  FixedMatrix<typename MulType<T, T1>::Type, row_num, new_col_num> operator*(const FixedMatrix<T1, col_num, new_col_num>& other) const {
    typedef typename MulType<T, T1>::Type ElementType;
    typedef FixedMatrix<ElementType, row_num, new_col_num> ResultMatrix;
    ResultMatrix result;
    FixedMatOp<T, row_num, col_num>::template Mul<T1, new_col_num>((T*) matrix_,
                                                                   (T1*) other.matrix_,
                                                                   (ElementType*) result.matrix_);
    return result;
  }

  template <class T1, int new_col_num>
  FixedMatrix<typename MulType<T, T1>::Type, row_num, new_col_num> operator*(const FixedMatrixWrapper<T1, col_num, new_col_num>& other) const {
    typedef typename MulType<T, T1>::Type ElementType;
    typedef FixedMatrix<ElementType, row_num, new_col_num> ResultMatrix;
    ResultMatrix result;
    FixedMatOp<T, row_num, col_num>::template Mul<T1, new_col_num>((T*) matrix_, (T1*) other.matrix_, (ElementType*) result.matrix_);
    return result;
  }


  template <class T1, bool wrapped>
  Vec<typename MulType<T, T1>::Type, row_num, false> operator*(const Vec<T1, col_num, wrapped>& vec) const {
    typedef typename MulType<T, T1>::Type ElementType;
    typedef Vec<ElementType, row_num, false> ResultVector;
    ResultVector result;
    FixedMatOp<T, row_num, col_num>::template Mul<T1>((T*) matrix_, vec(), result());
    return result;
  }

  template <class T1>
  Mat& operator*=(const T1& mul) {
//    static_assert(std::is_same<T, typename MulType<T, T1>::Type>::value,
//                 "Result type of multiplication must be the same as element type");
    FixedMatOp<T, row_num, col_num>::template Mul<T1>((T*) matrix_, mul, (T*) matrix_);
    return *this;
  }
  T matrix_[row_num][col_num];
};


template <class T, int row_num, int col_num>
inline FixedMatrix<typename MulType<T, T>::Type, row_num, col_num> operator*(const T& mul, const FixedMatrix<T, row_num, col_num>& m)
{
  typedef typename MulType<T, T>::Type ElementType;
  typedef FixedMatrix<ElementType, row_num, col_num> Mat;
  Mat result;
  FixedMatOp<T, row_num, col_num>::template MulElement<T>(m(), mul, result());
  return result;
}

// Disabled due to ambiguity
//template <class T, int row_num, int col_num>
//inline FixedMatrix<typename MulType<T, T>::Type, row_num, col_num> operator*(const FixedMatrix<T, row_num, col_num>& m, const T& mul)
//{
//  typedef typename MulType<T, T>::Type ElementType;
//  typedef FixedMatrix<ElementType, row_num, col_num> Mat;
//  Mat result;
//  FixedMatOp<T, row_num, col_num>::template MulElement<T>(m(), mul, result());
//  return result;
//}


template <class T, int row_num, int col_num, class T1, bool wrapped>
Vec<typename MulType<T, T1>::Type, col_num, false> operator*(const Vec<T1, row_num, wrapped>& vec,
                                                      const FixedMatrix<T, row_num, col_num>& matrix)
{
  typedef typename MulType<T, T1>::Type ElementType;
  typedef Vec<ElementType, col_num, wrapped> ResultVector;
  ResultVector result;
  FixedMatOp<T, row_num, col_num>::template MulLeft<T1>(vec(), matrix(), result());
  return result;
}

template <class T, int row_num, int col_num>
inline FixedMatrix<typename MulType<T, T>::Type, row_num, col_num> operator*(const T& mul, const FixedMatrixWrapper<T, row_num, col_num>& m)
{
  typedef typename MulType<T, T>::Type ElementType;
  typedef FixedMatrix<ElementType, row_num, col_num> Mat;
  Mat result;
  FixedMatOp<T, row_num, col_num>::template MulElement<T>(m(), mul, result());
  return result;
}

template <class T, int row_num, int col_num>
inline FixedMatrix<typename MulType<T, T>::Type, row_num, col_num> operator*(const FixedMatrixWrapper<T, row_num, col_num>& m, const T& mul)
{
  typedef typename MulType<T, T>::Type ElementType;
  typedef FixedMatrix<ElementType, row_num, col_num> Mat;
  Mat result;
  FixedMatOp<T, row_num, col_num>::template MulElement<T>(m(), mul, result());
  return result;
}

template <class T, int row_num, int col_num, class T1, bool wrapped>
Vec<typename MulType<T, T1>::Type, col_num, false> operator*(const Vec<T1, row_num, wrapped>& vec,
                                                      const FixedMatrixWrapper<T, row_num, col_num>& matrix)
{
  typedef typename MulType<T, T1>::Type ElementType;
  typedef Vec<ElementType, col_num, false> ResultVector;
  ResultVector result;
  FixedMatOp<T, row_num, col_num>::template MulLeft<T1>(vec(), matrix(), result());
  return result;
}

} // namespace;


#endif // FIXED_MATRIX_H_
