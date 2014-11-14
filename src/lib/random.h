#ifndef RANDOM_H_
#define RANDOM_H_
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

class Random
{
public:
  Random(void) {
    srand((unsigned int) time(NULL));
  }
  template <class T>
  T GetRandom(T max = 1.0, T min = 0.0) {
    T random = rand() * ((T) 1.0) / RAND_MAX;
    return random * (max - min) + min;
  }
  int GetRandom(int max, int min = 0) {
    int random = rand();
    return (random % (max - min)) + min;
  }
};
#endif
