#include "main_window.h"
#include <QDebug>
#include <QApplication>
#include <iterator>
#include <typeinfo>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include "global.h"
#include "vector_lib.h"
#include "print_macro.h"

void test()
{
}

int main(int argc, char *argv[])
{
#if 0
  test();
#else
  QApplication a(argc, argv);
  setlocale(LC_NUMERIC, "C");
  MainWindow w;
  w.show();
  return a.exec();
#endif
}

