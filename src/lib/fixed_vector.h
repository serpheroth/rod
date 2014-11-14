#ifndef FIXED_VECTOR_H_
#define FIXED_VECTOR_H_
#include <type_traits>
#include "utility_function.h"
#include "fixed_vector_utility.h"
#include "macro_constant.h"

#ifndef _MSC_VER
#define ENABLE_UNION
#endif
#undef ENABLE_UNION
//#define ENABLE_CXX_11
namespace dj
{

template <class ElementType, int d>
struct Storage {
  ElementType v[d];
  Storage() {}
  explicit Storage(const ElementType* v_in) { FixedVecOp<d, ElementType>::Copy(v_in, v); }
};

template <class ElementType>
struct Storage<ElementType, 2> {
#ifdef ENABLE_UNION
  union {
    struct { ElementType x, y; };
    ElementType v[2];
  };
  Storage(ElementType a, ElementType b) : x(a), y(b) {}
  explicit Storage(const ElementType* v_in) : x(v[0]), y(v[1]) {}
#else
  ElementType v[2];
# ifndef ENABLE_CXX_11
  Storage(ElementType a, ElementType b)  { v[0] = a; v[1] = b; }
  explicit Storage(const ElementType* v_in) { v[0] = v_in[0]; v[1] = v_in[1]; }
# else
  Storage(ElementType a, ElementType b) : v({a, b}) {}
  explicit Storage(const ElementType* v_in) : v({v_in[0], v_in[1]}) {}
# endif
#endif
  Storage() {}
};

template <class ElementType>
struct Storage<ElementType, 3> {
#ifdef ENABLE_UNION
  union {
    struct { ElementType x, y, z; };
    ElementType v[3];
  };
  Storage(ElementType a, ElementType b, ElementType c) : x(a), y(b), z(c) {}
  explicit Storage(const ElementType* v_in) : x(v[0]), y(v[1]), z(v[2]) {P(x, y, z);}
#else
  ElementType v[3];
# ifndef ENABLE_CXX_11
  Storage(ElementType a, ElementType b, ElementType c)  { v[0] = a; v[1] = b; v[2] = c; }
  explicit Storage(const ElementType* v_in) { v[0] = v_in[0]; v[1] = v_in[1]; v[2] = v_in[2]; }
# else
  Storage(ElementType a, ElementType b, ElementType c) : v({a, b, c}) {}
  explicit Storage(const ElementType* v_in) : v({v_in[0], v_in[1], v_in[2]}) {}
# endif
#endif
  Storage() {}
};

template <class ElementType>
struct Storage<ElementType, 4> {
#ifdef ENABLE_UNION
  union {
    struct { ElementType x, y, z, w; };
    ElementType v[4];
  };
  Storage(ElementType a, ElementType b, ElementType c, ElementType d) : x(a), y(b), z(c), w(d) {}
  explicit Storage(const ElementType* v_in) : x(v[0]), y(v[1]), z(v[2]), w(v[3]) {}
#else
  ElementType v[4];
#  ifndef ENABLE_CXX_11
  Storage(ElementType a, ElementType b, ElementType c, ElementType d)  { v[0] = a; v[1] = b; v[2] = c; v[3] = d; }
  explicit Storage(const ElementType* v_in) { v[0] = v_in[0]; v[1] = v_in[1]; v[2] = v_in[2]; v[2] = v_in[3]; }
#  else
  Storage(ElementType a, ElementType b, ElementType c, ElementType d) : v({a, b, c, d}) {}
  explicit Storage(const ElementType* v_in) : v({v_in[0], v_in[1], v_in[2], v_in[3]}) {}
#  endif
#endif
  Storage() {}
};

template <class ElementType, int d, bool wrapped = false>
struct StorageChooser {
  typedef Storage<ElementType, d> data_type;
};

template <class ElementType, int d>
struct StorageChooser<ElementType, d, true> {
  struct data_type {
    ElementType * v;
    data_type() {}
    explicit data_type(ElementType* v_in) : v(v_in) {}
  };
};

template <class T, int d, bool wrapped>
class Vec : public StorageChooser<T, d, wrapped>::data_type
{
public:
  typedef T* iterator;
  typedef const T* const_iterator;
  typedef typename StorageChooser<T, d, wrapped>::data_type data_type;
  //  static const int kDimension = d;
  enum {
    kDimension = d
  };
  enum {
    kWrapped = wrapped
  };
#ifndef _MSC_VER
  using data_type::v;
#endif
  explicit Vec(const T initial_value = T())
  {
    static_assert(kWrapped == false, "constructor not applicable to wrapped vector");
    FixedVecOp<d, T>::Fill(v, initial_value);
  }

  explicit Vec(T* initial_values) : data_type(initial_values) { }
  explicit Vec(const T* initial_values) : data_type(initial_values)
  {
    static_assert(wrapped == false, "constructor not applicable to wrapped vector");
  }
  Vec(T a, T b) : data_type(a, b) {}
  Vec(T a, T b, T c) : data_type(a, b, c) {}
  Vec(T a, T b, T c, T e) : data_type(a, b, c, e) {}

  /*
  Vec(std::initializer_list<T> list) {
    auto iter = list.begin();
    for (int i = 0; i < d; ++i, ++iter) {
      x[i] = *iter;
    }
  }*/

  template <bool input_wrapped>
  Vec(const Vec < T, d - 1, input_wrapped > & o)
  {
    static_assert(wrapped == false, "constructor not applicable to wrapped vector");
    FixedVecOp<d, T>::Copy(o.v, v);
    v[d - 1] = 0;
  }

  template <bool input_wrapped>
  Vec(const Vec < T, d - 1, input_wrapped > & o, T last_element)
  {
    static_assert(wrapped == false, "constructor not applicable to wrapped vector");
    FixedVecOp<d, T>::Copy(o.v, v);
    v[d - 1] = last_element;
  }

  template <bool input_wrapped>
  Vec(const Vec<T, d, input_wrapped>& other)
  {
    static_assert(wrapped == false, "constructor not applicable to wrapped vector");
    FixedVecOp<d, T>::Copy(other.v, v);
  }

  template <bool input_wrapped>
  const Vec<T, d, wrapped>& operator=(const Vec<T, d, input_wrapped>& other)
  {
    FixedVecOp<d, T>::Copy(other.v, v);
    return *this;
  }

  inline const T& operator[](int index) const
  {
    return v[index];
  }

  inline T& operator[](int index)
  {
    return v[index];
  }

  inline T& At(int index)
  {
    if (index < 0 || index >= kDimension) {
      std::cerr << "Vec:operator[] => index " << index << " is out of range [0, " << d - 1 << "]" << std::endl;
      return v[0];
    } else {
      return v[index];
    }
  }
  template <class VectorType>
  bool operator==(const VectorType& other) const
  {
    return FixedVecOp<d, T>::Equal(v, &other[0]);
  }

  template <class VectorType>
  bool operator!=(const VectorType& other) const
  {
    return !(FixedVecOp<d, T>::Equal(v, &other[0]));
  }

  template <class VectorType>
  bool operator<(const VectorType& other) const
  {
    return FixedVecOp<d, T>::LessThan(v, &other[0]);
  }

  template <class VectorType>
  bool operator>(const VectorType& other) const
  {
    return FixedVecOp<d, T>::LessThan(&other[0], v);
  }

  iterator begin() { return v; }
  const_iterator begin() const { return v; }
  iterator end() { return v + d; }
  const_iterator end() const { return v + d; }

  template <class VectorType>
  Vec<T, d, false> operator+(const VectorType& other) const
  {
    Vec<T, d, false> result;
    FixedVecOp<d, T>::Plus(v, &other[0], result.v);
    return result;
  }

  template <class VectorType>
  Vec<T, d, false> operator-(const VectorType& other) const
  {
    Vec<T, d, false> result;
    FixedVecOp<d, T>::Minus(v, &other[0], result.v);
    return result;
  }

  template <class VectorType>
  Vec<T, d, wrapped>& operator-=(const VectorType& other)
  {
    FixedVecOp<d, T>::Minus(v, &other[0], v);
    return *this;
  }

  template <class VectorType>
  Vec<T, d, wrapped>& operator+=(const VectorType& other)
  {
    FixedVecOp<d, T>::Plus(v, other.v, (T*) v);
    return *this;
  }

  Vec<T, d, wrapped>& operator/=(T mul)
  {
    FixedVecOp<d, T>::Scale(v, T(1) / mul);
    return *this;
  }

  Vec<T, d, wrapped>& operator*=(T mul)
  {
    FixedVecOp<d, T>::Scale(v, mul);
    return *this;
  }

  void Fill(T value)
  {
    std::fill(v, v + d, value);
  }

  T MagnitudeSquare(void) const
  {
    return FixedVecOp<d, T>::MagnitudeSquare(v);
  }

  T Magnitude(void) const
  {
    return FixedVecOp<d, T>::Magnitude(v);
  }

  T Normalize(void)
  {
    return FixedVecOp<d, T>::Normalize(v);
  }

  inline T* operator()() { return v; }
  inline const T* operator()() const { return v;}
  inline T* Data() { return v; }
  inline const T* Data() const { return v;}

  template <bool input_wrapped>
  T DistanceTo(const Vec<T, d, input_wrapped>& to) const
  {
    return FixedVecOp<d, T>::Distance(v, to.v);
  }

  template <bool input_wrapped>
  T DistanceSquareTo(const Vec<T, d, input_wrapped>& to) const
  {
    return FixedVecOp<d, T>::DistanceSquare(v, to.v);
  }

  T Max(void) const
  {
    return FixedVecOp<d, T>::Max(v);
  }

  T Min(void) const
  {
    return FixedVecOp<d, T>::Min(v);
  }
};

template <class T1, int d, class T2, bool arg1_wrapped, bool arg2_wrapped>
typename MulType<T1, T2>::Type operator*(const Vec<T1, d, arg1_wrapped>& v, const Vec<T2, d, arg2_wrapped>& o)
{
  return FixedVecOp<d, T1>::Dot(&v[0], o.v);
}

template <class T, int d, class T1, bool wrapped>
inline Vec<T, d, wrapped> operator*(const T1 mul, const Vec<T, d, wrapped>& x)
{
  Vec<T, d, false> result = x;
  FixedVecOp<d, T>::Scale(result.v, mul);
  return result;
}

template <class T, int d, class T1, bool wrapped>
inline Vec<T, d, false> operator*(const Vec<T, d, wrapped>& x, const T1 mul)
{
  Vec<T, d, false> result = x;
  FixedVecOp<d, T>::Scale(result.v, mul);
  return result;
}

template <class T, int d, class T1, bool wrapped>
inline Vec<T, d, false> operator/(const Vec<T, d, wrapped>& x, const T1 mul)
{
  return operator*(T(1) / mul, x);
}

template <class T, int d, bool wrapped>
std::ostream& operator<<(std::ostream& out, const Vec<T, d, wrapped>& vec)
{
  out << "[" << vec.v[0];
  for (int i = 1; i < d; i++) {
    out << ", " << vec.v[i];
  }
  out << "]";
  return out;
}

typedef Vec<unsigned int, 2, false> Vec2ui;
typedef Vec<unsigned int, 3, false> Vec3ui;
typedef Vec<unsigned int, 4, false> Vec4ui;
typedef Vec<int, 2, false> Vec2i;
typedef Vec<int, 3, false> Vec3i;
typedef Vec<int, 4, false> Vec4i;
typedef Vec<float, 2, false> Vec2f;
typedef Vec<float, 3, false> Vec3f;
typedef Vec<float, 4, false> Vec4f;
typedef Vec<double, 2, false> Vec2d;
typedef Vec<double, 3, false> Vec3d;
typedef Vec<double, 4, false> Vec4d;
}
#endif // FIXED_VECTOR_H_
