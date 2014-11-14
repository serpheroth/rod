#include <iostream>
#include <algorithm>
#include "blas_wrapper.h"
#include "conjugate_gradient_solver.h"
#include "print_macro.h"
const int ConjugateGradientSolver::kDefaultMaxIteration = 200;
const Float ConjugateGradientSolver::kDefaultResidualThreadshold = 1e-16;
ConjugateGradientSolver::ConjugateGradientSolver(int size) : size_(size)
{
#ifdef MKL
  mkl_set_num_threads(4);
#endif
  normalizer_ = 1.0 / size;
  residual_ = buffer_manager.Malloc<Float>(size_);
  search_direction_ = buffer_manager.Malloc<Float>(size_);
  temp_ = buffer_manager.Malloc<Float>(size_);
}

inline void ConjugateGradientSolver::ForEach(Float* v1, Float* v2, Float* out, Transformer transform)
{
  for (int i = 0; i < size_; i++) {
    if (filter_(i)) {
      out[i] = transform(v1[i], v2[i]);
    }
  }
}

inline Float ConjugateGradientSolver::Dot(Float* v1, Float* v2)
{
  Float dot = 0;
  for (int i = 0; i < size_; i++) {
    if (filter_(i)) {
      dot += v1[i] * v2[i];
    }
  }
  return dot;
}

ConjugateGradientSolver::Info ConjugateGradientSolver::Solve(Float* rhs,
                                                             Float* result,
                                                             Matrix A,
                                                             Filter filter,
                                                             int max_iteration,
                                                             Float residual_threshold)
{
  Float current_residual = 0;
  max_iteration = (max_iteration == 0) ? kDefaultMaxIteration : max_iteration;
  filter_ = filter;
  memset(temp_, 0, sizeof(Float) * size_);
  memset(result, 0, sizeof(Float) * size_); // result = 0
  memcpy(residual_, rhs, sizeof(Float) * size_); // r = rhs - A*0 = rhs
  memcpy(search_direction_, residual_, sizeof(Float) * size_);  // p = rhs
  Float denumerator;
  denumerator = Dot(residual_, residual_);
  float residual_2_norm = denumerator;
  int iteration = 0 ;
  if (denumerator < -0.0000001 || denumerator > 0.00000001) {
    for (iteration = 1; iteration <= max_iteration; iteration++) {
      A(search_direction_, temp_) ; // temp = A * p
            Float alpha = Dot(search_direction_, temp_);
      alpha = denumerator / alpha ;  // alpha = r*r/p*A*p
      ForEach(search_direction_, result, result, [&](Float a, Float b) { return alpha * a + b; } ); // result += alpha * search_direction_;
      ForEach(temp_, residual_, residual_, [&](Float a, Float b) { return -alpha * a + b; } ); // residual_ -= alpha * A * search_direction_;
      residual_2_norm = Dot(residual_, residual_);
      if (residual_2_norm != residual_2_norm) {
        std::cerr << CURRENT_LINE << " => " << iteration << " iterations with residual " << residual_2_norm << std::endl;
        getchar();
      }
      current_residual = residual_2_norm * normalizer_;
      if (current_residual < residual_threshold) {
        //std::cout << CURRENT_LINE << " => " << iteration << " iterations with residual " << residual_2_norm <<std::endl;
        break ;
      }
      Float beta = residual_2_norm / denumerator;
      ForEach(residual_, search_direction_, search_direction_, [&](Float a, Float b) { return a + b * beta;} );  // p = r+beta*p
      denumerator = residual_2_norm;
    }
    if (iteration > max_iteration) {
      std::cerr << CURRENT_LINE << " => warning: maximum " << max_iteration << " iteration reached with residual " << residual_2_norm << std::endl;
    }
  }
  return Info(iteration, current_residual) ;
}

ConjugateGradientSolver::Info ConjugateGradientSolver::Solve(Float* rhs,
                                                             Float* result,
                                                             Matrix A,
                                                             int max_iteration,
                                                             Float residual_threshold)
{
  Float current_residual = 0;
  max_iteration = (max_iteration == 0) ? kDefaultMaxIteration : max_iteration;
//  memset(result, 0, sizeof(Float) * size_); // result = 0
  A(result, temp_);
  memcpy(residual_, rhs, sizeof(Float) * size_);  // p = rhs
  blas::axpy(size_, -1.0f, temp_, residual_);
  memcpy(search_direction_, residual_, sizeof(Float) * size_);  // p = rhs
  Float denumerator;
  blas::dot(size_, residual_, residual_, denumerator);
  float residual_2_norm = denumerator;
  int iteration = 0 ;
  ASSERT(denumerator == denumerator);
  if (denumerator * normalizer_ > residual_threshold) {
    for (iteration = 1; iteration <= max_iteration; iteration++) {
      A(search_direction_, temp_) ; // temp = A * p
      Float alpha;
      blas::dot(size_, search_direction_, temp_, alpha);
      alpha = denumerator / alpha ;  // alpha = r*r/p*A*p
      blas::axpy(size_, alpha, search_direction_, result); // x = x+alpha*p
      blas::axpy(size_, -alpha, temp_, residual_);  // r = r-alpha*A*p
      blas::dot(size_, residual_, residual_, residual_2_norm) ;  // residual = r*r
      if (residual_2_norm != residual_2_norm) {
        std::cerr << CURRENT_LINE << " => " << iteration << " iterations with residual " << residual_2_norm << std::endl;
        getchar();
      }
      current_residual = residual_2_norm * normalizer_;
      if (current_residual < residual_threshold) {
        std::cout << CURRENT_LINE << " => " << iteration << " iterations with residual " << current_residual <<std::endl;
        break ;
      }
      Float beta = residual_2_norm / denumerator;
      blas::scal(size_, beta, search_direction_) ; // p = r+beta*p
      blas::axpy(size_, 1.0 , residual_, search_direction_);
      denumerator = residual_2_norm;
    }
    if (iteration > max_iteration) {
      std::cerr << CURRENT_LINE << " => warning: maximum " << max_iteration << " iteration reached with residual " << residual_2_norm << std::endl;
    }
  }
  return Info(iteration, current_residual) ;
}
