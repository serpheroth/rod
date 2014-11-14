#ifndef BLAS_WRAPPER_H_
#define BLAS_WRAPPER_H_
//#define MKL
#ifdef MKL
#include <mkl_types.h>
#include <mkl_cblas.h>
#include <mkl_service.h>
#else
#include <gsl/gsl_blas.h>
#endif
namespace blas
{
// vec2 += coefficient*vec1
inline void axpy(int size,float coefficient,float* vec1, float* vec2)
{
  cblas_saxpy(size, coefficient,vec1, 1, vec2, 1); // x = x+alpha*p
}

// vec2 += coefficient*vec1
inline void axpy(int size,double coefficient,double* vec1, double* vec2)
{
  cblas_daxpy(size, coefficient,vec1, 1, vec2, 1); // x = x+alpha*p
}

// dot_result = vec1*vec2
inline void dot(int size,float* vec1, float* vec2,float& dot_result)
{
  dot_result = cblas_sdot(size, &vec1[0], 1, &vec2[0], 1);
}
inline void dot(int size,double* vec1, double* vec2,double& dot_result)
{
  dot_result = cblas_ddot(size, &vec1[0], 1, &vec2[0], 1);
}
inline double dotd(int size,double* vec1, double* vec2)
{
  return cblas_ddot(size, &vec1[0], 1, &vec2[0], 1);
}
inline float dotf(int size,float* vec1, float* vec2)
{
  return cblas_sdot(size, &vec1[0], 1, &vec2[0], 1);
}

// vec = coefficient*vec
inline void scal(int size,float coefficient,float* vec)
{
  cblas_sscal(size, coefficient, vec, 1);
}
inline void scal(int size,double coefficient,double* vec)
{
  cblas_dscal(size, coefficient, vec, 1);
}

}
#endif
