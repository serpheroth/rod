#include "rainbow_color.h"
//#include "print_macro.h"
#include "vector_lib.h"
/*
float spectrum[][3] = {
{1.0f,0.5f,0.5f},{1.0f,0.75f,0.5f},{1.0f,1.0f,0.5f },{0.75f,1.0f,0.5f},
{0.5f,1.0f,0.5f},{0.5f,1.0f,0.75f},{0.5f,1.0f,1.0f },{0.5f,0.75f,1.0f},
{0.5f,0.5f,1.0f},{0.75f,0.5f,1.0f},{1.0f,0.5f,1.0f },{1.0f,0.5f,0.75f}
};

float spectrum[][3] = {
0,0,1, //blue
0.56078431372f, 0, 255, //violet
0.29411764705, 0, 0.50980392156, //indigo
0,1,0, // green
1,1,1, // yellow
1,0.49803921568,0, // orange
1,0,0 // red
};
*/
float spectrum[][3] = {
  {0,0,1}, //blue
//0.56078431372f, 0, 255, //violet
//0.29411764705, 0, 0.50980392156, //indigo
  {0,1,0}, // green
  {1,1,0}, // yellow
//1,0.49803921568,0, // orange
  {1,0,0} // red
};
const dj::Vec3f kBlue(0.00f, 0.00f, 1.00f);
const dj::Vec3f kViolet(0.56f, 0.00f, 1.00f);
const dj::Vec3f kIndigo(0.29f, 0.00f, 0.51f);
const dj::Vec3f kGreen(0.00f, 1.00f, 0.00f);
const dj::Vec3f kYellow(1.00f, 1.00f, 0.00f);
const dj::Vec3f kOrage(1.00f, 0.50f, 0.00f);
const dj::Vec3f kRed(1.00f, 0.00f, 0.00f);
const dj::Vec3f kGrey(0.59f, 0.59f, 0.59f);
const dj::Vec3f kBlack(0.00f, 0.00f, 0.00f);
const dj::Vec3f kWhite(1.00f, 1.00f, 1.00f);
const dj::Vec3f kChocolate(0.188f, 0.0980f, 0.067f);


float GetColor(float value, float* color)
{
  if ( value != value ) {
    color[0]=color[1]=color[2] = 1.0 ;
    std::cerr << "getColor(): input value is nan" << std::endl ;
    return 1;
  }
  if ( value >=1 ) {
    color[0]=color[1]=color[2] = 1.0 ;
    return 1;
  }
  if ( value <=0 ) {
    color[0]=color[1]= 0.0 ;
    color[2] =1.0 ;
    return 0;
  }

  int left = int(value*3) ;
  float ratio = value*3 - left ;
  if ( left >= 3 || left<0 ) {
    std::cerr<<"getColor(): invalid index value: "<<value<<std::endl;
  }
  color[0] = (1-ratio)*spectrum[left][0] + ratio*spectrum[left+1][0] ;
  color[1] = (1-ratio)*spectrum[left][1] + ratio*spectrum[left+1][1] ;
  color[2] = (1-ratio)*spectrum[left][2] + ratio*spectrum[left+1][2] ;
  return value ;
}

float GetColor(float min, float max, float value, float* color)
{
  double gap = max - min ;
  if ( gap<0.0001 ) {
    color[0]=spectrum[0][0];
    color[1]=spectrum[0][1];
    color[2]=spectrum[0][2];
    return 0;
  }
  value = (value - min ) / (max-min) ;
  if ( value != value ) {
    color[0]=color[1]=color[2] = 1.0 ;
    std::cerr << "getColor(): input value is nan" << std::endl ;
    return 1;
  }

  if ( value >=1 ) {
    color[0]=color[1]=color[2] = 0.0 ;
    color[0] = 1.0 ;
    return 1;
  }
  if ( value <=0 ) {
    color[0]=color[1]= 0.0 ;
    color[2] =1.0 ;
    return 0;
  }
  int left = int(value*3) ;
  float ratio = value*3 - left ;
  if ( left >= 3 || left<0 ) {
    std::cerr<<"getColor(): invalid index value: "<<value<<std::endl;
  }
  color[0] = (1-ratio)*spectrum[left][0] + ratio*spectrum[left+1][0] ;
  color[1] = (1-ratio)*spectrum[left][1] + ratio*spectrum[left+1][1] ;
  color[2] = (1-ratio)*spectrum[left][2] + ratio*spectrum[left+1][2] ;
  return value ;
}
