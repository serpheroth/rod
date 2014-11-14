#ifndef WU_PRINT_MACRO_H_
#define WU_PRINT_MACRO_H_
#include <iostream>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include "macro_constant.h"

//#define ENABLE_COLORED_TEXT_OUTPUT
//#define TWO_COLOR_OUT_PUT
//********************************************************************************
//--------------------------------------------------------------------------------
//	 Macros for printing out variable for the purpose of debuging
//--------------------------------------------------------------------------------
//********************************************************************************

#if defined(ENABLE_COLORED_TEXT_OUTPUT) && !defined(TWO_COLOR_OUT_PUT)
#define BLACK     "\033[1;30m"
#define RED       "\033[1;31m"
#define GREEN     "\033[1;32m"
#define YELLOW    "\033[1;33m"
#define BLUE      "\033[1;34m"
#define MAGENTA   "\033[1;35m"
#define CYAN      "\033[1;36m"
#define WHITE     "\033[1;37m"
#define NORMAL    "\033[0m"
#else
#define BLACK     ""
#define RED       ""
#define GREEN     ""
#define YELLOW    ""
#define BLUE      ""
#define MAGENTA   ""
#define CYAN      ""
#define WHITE     ""
#define NORMAL    ""
#endif // ENABLE_COLORED_TEXT_OUTPUT

#define COLOR(x,y) (y) << (x) << NORMAL
template <class T>
inline std::string SetColor(T value, const char* color = RED)
{
  std::stringstream colored_string;
  colored_string << color << value << NORMAL;
  return colored_string.str();
}

#define VALUE_COLOR(x) GREEN << (x) << NORMAL
template <class T>
inline void PrintToken(const char* token, T value) {
  std::cout << token << " = ";
  std::cout.flush();
  std::cerr << value;
  std::cerr.flush();
}

#define MAX_FILE_NAME_LENGTH 25
#define MAX_LINE_NUM_LENGTH 4
#define MAX_FUNCTION_NAME_LENGTH 30
//#define CURRENT_LINE            (std::string(__FILE__) + std::string(":") + std::string(ToString(__LINE__)) + std::string(":") +
#ifdef TWO_COLOR_OUT_PUT
#define PRINT_CURRENT_LINE      std::cerr << MAGENTA << std::setw(MAX_FILE_NAME_LENGTH) << std::right << __FILE__ << " " << std::setw(MAX_LINE_NUM_LENGTH) <<   __LINE__  << " " << std::setw(MAX_FUNCTION_NAME_LENGTH) << std::left << FUNCTION_NAME  << NORMAL << " => ";
#else
#define PRINT_CURRENT_LINE      std::cout << MAGENTA << std::setw(MAX_FILE_NAME_LENGTH) << std::right << __FILE__ << " " << std::setw(MAX_LINE_NUM_LENGTH) <<   __LINE__  << " " << std::setw(MAX_FUNCTION_NAME_LENGTH) << std::left << FUNCTION_NAME  << NORMAL << " => ";
#endif

template <class T>
inline std::string ToString(const char* token, T value)
{
  std::stringstream str;
  str << token << " = " << VALUE_COLOR(value);
  return str.str();
}
//--------------------------------------------------------------------------------
//	print out list the variable names and value seperated with tab
//--------------------------------------------------------------------------------
#ifndef TWO_COLOR_OUT_PUT
#define PL1(x)		  {\
 std::string STRING = ToString(STRINGIZE_TOKEN(x), x); \
  PRINT_CURRENT_LINE;\
  std::cout << STRING;\
}
#define PL2(x,y)	{\
  std::string STRING1 = ToString(STRINGIZE_TOKEN(x), x); \
  std::string STRING2 = ToString(STRINGIZE_TOKEN(y), y); \
  PRINT_CURRENT_LINE;  \
  std::cout << STRING1 << "\t" << STRING2; \
}
#define PL3(x,y,z) {\
  std::string STRING1 = ToString(STRINGIZE_TOKEN(x), x); \
  std::string STRING2 = ToString(STRINGIZE_TOKEN(y), y); \
  std::string STRING3 = ToString(STRINGIZE_TOKEN(z), z); \
  PRINT_CURRENT_LINE; \
  std::cout << STRING1 << "\t" << STRING2 << "\t" << STRING3; \
}

#define PL4(x, y, z, a) {\
  std::string STRING1 = ToString(STRINGIZE_TOKEN(x), x); \
  std::string STRING2 = ToString(STRINGIZE_TOKEN(y), y); \
  std::string STRING3 = ToString(STRINGIZE_TOKEN(z), z); \
  std::string STRING4 = ToString(STRINGIZE_TOKEN(a), a); \
  PRINT_CURRENT_LINE; \
  std::cout << STRING1 << "\t" << STRING2 << "\t" << STRING3 << "\t" << STRING4; \
}
#else
#define PL1(x)		  {
std::string STRING = ToString(STRINGIZE_TOKEN(x), x); \
PRINT_CURRENT_LINE;\
PrintToken(STRINGIZE_TOKEN(x), x);\
}

#define PL2(x,y)	{\
	PRINT_CURRENT_LINE;\
	PrintToken(STRINGIZE_TOKEN(x), x); \
	std::cout << "\t"; \
	PrintToken(STRINGIZE_TOKEN(y), y); \
}
#define PL3(x,y,z) {\
	PRINT_CURRENT_LINE;\
	PrintToken(STRINGIZE_TOKEN(x), x); \
	std::cout << "\t"; \
	PrintToken(STRINGIZE_TOKEN(y), y); \
	std::cout << "\t"; \
	PrintToken(STRINGIZE_TOKEN(z), z); \
}

#define PL4(x, y, z, a) {\
	PRINT_CURRENT_LINE;\
	PrintToken(STRINGIZE_TOKEN(x), x); \
	std::cout << "\t"; \
	PrintToken(STRINGIZE_TOKEN(y), y); \
	std::cout << "\t"; \
	PrintToken(STRINGIZE_TOKEN(z), z); \
	std::cout << "\t"; \
	PrintToken(STRINGIZE_TOKEN(a), a); \
}
#endif

//#define PL_STRING_MACRO_CHOOSER(...) GET_5TH_ARG(__VA_ARGS__, PL4, PL3, PL2, PL1, )
#define PL_CHOOSE_HELPER4(count) PL##count
#define PL_CHOOSE_HELPER3(count) PL_CHOOSE_HELPER4(count)
#define PL_CHOOSE_HELPER2(count) PL_CHOOSE_HELPER3(count)
#define PL_CHOOSE_HELPER1(count) PL_CHOOSE_HELPER2(count)
#define PL_CHOOSE_HELPER(count)  PL_CHOOSE_HELPER1(count)
//#define PL_STRING_MACRO_CHOOSER(...) GET_5TH_ARG(__VA_ARGS__, PL4, PL3, PL2, PL1, )
//#define PL(...) PL_STRING_MACRO_CHOOSER(__VA_ARGS__)(__VA_ARGS__)
#define PL(...) \
  MY_GLUE(PL_CHOOSE_HELPER(COUNT_ARGS(__VA_ARGS__)), (__VA_ARGS__))
//--------------------------------------------------------------------------------
//	print out list the variable names and value seperated with tab and end with new line
//	e.g. int a = 1, b = 2;
//      P(a,b) will print out "a=1	b=2"
//--------------------------------------------------------------------------------
#define P1(x)		  {\
	PL1(x); \
	std::cout << std::endl; \
}

#define P2(x,y)		{\
	PL2(x, y); \
	std::cout << std::endl; \
}

#define P3(x,y,z) {\
  PL3(x, y, z); \
  std::cout << std::endl; \
}

#define P4(x,y,z,a) {\
  PL4(x, y, z, a); \
  std::cout << std::endl; \
}

//#define P_STRING_MACRO_CHOOSER(...) GET_5TH_ARG(__VA_ARGS__, P4, P3, P2, P1, )
#define P_CHOOSE_HELPER4(count) P##count
#define P_CHOOSE_HELPER3(count) P_CHOOSE_HELPER4(count)
#define P_CHOOSE_HELPER2(count) P_CHOOSE_HELPER3(count)
#define P_CHOOSE_HELPER1(count) P_CHOOSE_HELPER2(count)
#define P_CHOOSE_HELPER(count)  P_CHOOSE_HELPER1(count)
#define P(...) \
  MY_GLUE(P_CHOOSE_HELPER(COUNT_ARGS(__VA_ARGS__)), (__VA_ARGS__))

//#define P(...) P_STRING_MACRO_CHOOSER(__VA_ARGS__)(__VA_ARGS__)
//--------------------------------------------------------------------------------
//	print out the values separated by tab
//--------------------------------------------------------------------------------
#define LL1(x)		  std::cout<< RED << (x) << NORMAL <<"\t"
#define LL2(x,y)		LL1(x); LL1(y)
#define LL3(x,y,z)	LL1(x); LL1(y); LL1(z)
#define LL4(x,y,z, a)	LL1(x); LL1(y); LL1(z); LL1(a)

#define LL_CHOOSE_HELPER4(count) LL##count
#define LL_CHOOSE_HELPER3(count) LL_CHOOSE_HELPER4(count)
#define LL_CHOOSE_HELPER2(count) LL_CHOOSE_HELPER3(count)
#define LL_CHOOSE_HELPER1(count) LL_CHOOSE_HELPER2(count)
#define LL_CHOOSE_HELPER(count)  LL_CHOOSE_HELPER1(count)
#define LL(...) \
  MY_GLUE(LL_CHOOSE_HELPER(COUNT_ARGS(__VA_ARGS__)), (__VA_ARGS__))
//#define LL_STRING_MACRO_CHOOSER(...) GET_5TH_ARG(__VA_ARGS__, LL4, LL3, LL2, LL1, )
//#define LL(...) LL_STRING_MACRO_CHOOSER(__VA_ARGS__)(__VA_ARGS__)
//--------------------------------------------------------------------------------
//	print out the values separated by tab and start a new line
//--------------------------------------------------------------------------------
#define L1(x)		      PRINT_CURRENT_LINE; LL1(x);     NEW_LINE
#define L2(x,y)	  	  PRINT_CURRENT_LINE; LL2(x,y);   NEW_LINE
#define L3(x,y,z)		  PRINT_CURRENT_LINE; LL3(x,y,z); NEW_LINE
#define L4(x,y,z,a)		  PRINT_CURRENT_LINE; LL4(x,y,z,a); NEW_LINE
//#define L_STRING_MACRO_CHOOSER(...) GET_5TH_ARG(__VA_ARGS__, L4, L3, L2, L1, )
//#define L(...) L_STRING_MACRO_CHOOSER(ME_PASS_VA(__VA_ARGS__))(ME_PASS_VA(__VA_ARGS__))
//#define L(...) L1(STRINGIZE_TOKEN(__VA_ARGS__))
//#define L(...) L1(ME_PASS_VA(__VA_ARGS__))

#define L_CHOOSE_HELPER4(count) L##count
#define L_CHOOSE_HELPER3(count) L_CHOOSE_HELPER4(count)
#define L_CHOOSE_HELPER2(count) L_CHOOSE_HELPER3(count)
#define L_CHOOSE_HELPER1(count) L_CHOOSE_HELPER2(count)
#define L_CHOOSE_HELPER(count)  L_CHOOSE_HELPER1(count)
//#define PL_STRING_MACRO_CHOOSER(...) GET_5TH_ARG(__VA_ARGS__, PL4, PL3, PL2, PL1, )
//#define PL(...) PL_STRING_MACRO_CHOOSER(__VA_ARGS__)(__VA_ARGS__)
//#define L(...) L_CHOOSE_HELPER(COUNT_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define L(...) \
  MY_GLUE(L_CHOOSE_HELPER(COUNT_ARGS(__VA_ARGS__)), (__VA_ARGS__))

//--------------------------------------------------------------------------------
//	print out  a sequence of symbole to seperate the output
inline void B(void)
{
  for (int i = 0; i < 50; i++) std::cout << "-";
  std::cout << std::endl;
}

inline void B(const char x)
{
  for (int i = 0; i < 50; i++) std::cout << x;
  std::cout << std::endl;
}
inline void B(const char* x)
{
  for (int i = 0; i < 30; i++) std::cout << x;
  std::cout << std::endl;
}
#define KK	{std::cout << VALUE_COLOR(CURRENT_LINE);for(int i=0;i<50;i++) std::cout<< RED << "-"; std::cout<< NORMAL << std::endl;}
#define OL(x,y)  (std::cout << (x) <<"="<<y<<"\t")
#define O(x,y) (std::cout << (x) <<"="<<y<<std::endl)
//#define B(x) (std::cout << "\n====== " << x << " ======\n" << std::endl)
#define  E(x) (std::cerr << x << " <= "<< __FILE__ << ":" << FUNCTION_NAME << ":" << __LINE__ << std::endl )
// end of printing macro
//********************************************************************************

#endif // WU_PRINT_MACRO_H_
