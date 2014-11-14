#ifndef RAINBOW_COLOR_H_
#define RAINBOW_COLOR_H_
#include "vector_lib.h"
extern const dj::Vec3f kBlue;
extern const dj::Vec3f kViolet;
extern const dj::Vec3f kIndigo;
extern const dj::Vec3f kGreen;
extern const dj::Vec3f kYellow;
extern const dj::Vec3f kOrage;
extern const dj::Vec3f kRed;
extern const dj::Vec3f kGrey;
extern const dj::Vec3f kBlack;
extern const dj::Vec3f kWhite;
extern const dj::Vec3f kChocolate;

float GetColor(float value, float* color);
float GetColor(float min, float max, float value, float* color);
#endif // RAINBOW_COLOR_H
