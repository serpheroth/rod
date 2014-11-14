#ifndef GLOBAL_H_
#define GLOBAL_H_
#include <vector>
#include <string>
//#define OPEN_MP
#ifdef OPEN_MP
#  include <omp.h>
#  ifndef _MSC_VER
#  define OMP_FOR _Pragma("omp parallel for")
#  else
#  define OMP_FOR __pragma(omp parallel for)
#  endif // vs compiler
#else // no omp
#  define OMP_FOR
#endif // #ifndef OPEN_MP

#ifdef _MSC_VER
#undef _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
//#define _CRT_SECURE_NO_DEPRECATE
#endif // _MSC_VER

#if defined(_WIN32) || defined(_WIN64)
#define TOKEN_TO_STRING1(t) # t
#define STRINGIZE_TOKEN1(t) TOKEN_TO_STRING1(t)
#define DATA_DIRECTORY STRINGIZE_TOKEN1(CURRENT_DIRECTORY)
#endif
class OpenGLQt;
typedef double real;

namespace global
{
extern OpenGLQt* gl;
extern int simulate;
extern int pause_per_frame;
extern const char kPWD[];
extern const char kDataDir[];
extern std::string data_directory;
extern real time_step;
extern real* gravity;
extern int simulation_step_per_idle;
}

namespace rod {
extern real* bending_twising_stiffness;
extern int rod_v_num;
extern real rod_length;
extern int pbd_iteration;
}

#endif // GLOBAL_H_
