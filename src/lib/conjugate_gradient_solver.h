#ifndef CONJUGATE_GRADIENT_SOLVER_H_
#define CONJUGATE_GRADIENT_SOLVER_H_
#include <iostream>
#include <functional>
#include "blas_wrapper.h"
#include "buffer_manager.h"
#include "print_macro.h" // CURRENT_LINE
//#include "macro.h" // for typedef float Float
typedef float Float;
class ConjugateGradientSolver
{
 public:
  typedef std::function<Float (Float, Float)> Transformer;
  typedef std::function<bool (int)> Filter;
  typedef std::function<void (Float*, Float*)> Matrix;
  typedef std::pair<int, Float> Info; // <iteration to solve, residual>
 private:
  Float *residual_, *search_direction_, *temp_ ;
  Float normalizer_;
  Filter filter_;
  int size_ ;
  BufferManager buffer_manager;
  // Disallow assignment of one solver to another
  ConjugateGradientSolver();
  ConjugateGradientSolver(ConjugateGradientSolver&);
  const ConjugateGradientSolver& operator=(ConjugateGradientSolver&);
  inline void ForEach(Float* v1, Float* v2, Float* out, Transformer transform);
  inline Float Dot(Float* v1, Float* v2);
 public:
  static const int kDefaultMaxIteration;
  static const Float kDefaultResidualThreadshold;

  ConjugateGradientSolver(int size);

  ~ConjugateGradientSolver() {}

  void Resize(int size) {
    size_ = size;
  }

  Float set_normalizer(int n) {
    normalizer_ = 1.0 / n;
    return normalizer_;
  }

  // Compute residual = rhs - A * x
  inline void ComputeResidual(Float* x, Float* rhs, Matrix A, Float* residual) {
    A(x, residual);
    for (int i = 0; i < size_; i++) {
      residual[i] = rhs[i] - residual[i];
    }
  }

  Float ComputeResidual2Norm(Float* x, Float* rhs, Matrix A) {
    ComputeResidual(x, rhs, A, temp_);
    Float residual_2_norm;
    blas::dot(size_, temp_, temp_, residual_2_norm);
    return residual_2_norm;
  }

  //	A * result = rhs
  Info Solve(Float* rhs,
             Float* result,
             Matrix A,
             Filter filter,
             int max_iteration = kDefaultMaxIteration,
             Float residual_threshold = kDefaultResidualThreadshold);
  //	A * result = rhs
  Info Solve(Float* rhs,
             Float* result,
             Matrix A,
             int max_iteration = kDefaultMaxIteration,
             Float residual_threshold = kDefaultResidualThreadshold);
};
#endif // ConjugateGradientSolver
