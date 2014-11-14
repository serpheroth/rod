#ifndef MACRO_CONSTANT_H
#define MACRO_CONSTANT_H
#include <sstream>
#ifndef EPSILON
#define EPSILON 1e-15
#endif
#ifndef PI
#define PI 3.1415926f
#endif

#if defined(_WIN32) ||  defined(_WIN64)
#define FUNCTION_NAME __FUNCTION__
#define FORCE_INLINE __forceinline
#else
#define FUNCTION_NAME __func__
#define FORCE_INLINE __attribute__((always_inline))
#endif

#define UNUSED(x) ((void) x)

#define GET_5TH_ARG(arg1, arg2, arg3, arg4, arg5, ...) arg5
#define COUNT_ARGS_IMPL4(_1, _2, _3, _4, count, ...) count
#define COUNT_ARGS_IMPL(args) COUNT_ARGS_IMPL4 args
#define COUNT_ARGS(...) COUNT_ARGS_IMPL((__VA_ARGS__, 4, 3, 2, 1, 0))
#define MY_GLUE(x, y) x y

#define MY_DEBUG
#ifdef MY_DEBUG
#define MY_ASSERT1(x) \
  {\
     if (!(x)) { \
      std::cerr << CURRENT_LINE << " => " << STRINGIZE_TOKEN(x) << " failed!!!" << std::endl; \
      exit(-1); \
    } \
  }
#define MY_ASSERT2(x, y) \
  {\
     if (!(x)) { \
      y; \
      std::cerr << CURRENT_LINE << " => " << STRINGIZE_TOKEN(x) << " failed!!!" << std::endl; \
      exit(-1); \
    } \
  }
#else
#define MY_ASSERT1(x) {}
#define MY_ASSERT2(x, y) {}
#endif

#define MY_ASSERT_CHOOSE_HELPER2(count) MY_ASSERT##count
#define MY_ASSERT_CHOOSE_HELPER1(count) MY_ASSERT_CHOOSE_HELPER2(count)
#define MY_ASSERT_CHOOSE_HELPER(count)  MY_ASSERT_CHOOSE_HELPER1(count)
#define ASSERT(...) \
  MY_GLUE(MY_ASSERT_CHOOSE_HELPER(COUNT_ARGS(__VA_ARGS__)), (__VA_ARGS__))

#define EXECUTE_TIMES(a) \
{\
  static int count = 0; \
  if (count >= (a)) return; \
  else count++; \
}
#define TOKEN_TO_STRING(t) # t
#define STRINGIZE_TOKEN(t) TOKEN_TO_STRING(t)
#define NL	std::cout<<std::endl
#define NEW_LINE std::cout<<std::endl;
#ifndef OUT
#define OUT std::cout
#endif
#define ERR std::cerr
#define EL std::endl

template <class T>
inline std::string ToString(T value)
{
  std::stringstream ret_string;
  ret_string << value;
  return ret_string.str();
}


#define CURRENT_LINE (std::string(__FILE__) + std::string(":") + std::string(ToString(__LINE__)) + std::string(":") + std::string(FUNCTION_NAME))
#endif // MACRO_CONSTANT_H
