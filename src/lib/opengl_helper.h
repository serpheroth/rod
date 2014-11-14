#ifndef OPENGL_HELPER_H_
#define OPENGL_HELPER_H_
#if defined(_WIN32) || defined(_WIN64) || defined(WIN32) || defined(WIN64)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <cstddef>
#include <stdlib.h>
#ifdef __APPLE__
#include <GL/glew.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include <cmath>



int WritePPM(const char *output_file_name, unsigned char *image, int resolution_x, int resolution_y);
/**
 * @brief GetPointDepth given a point in its local space, return its depth in NDC space
 */
double GetPointDepth(double *pos);
float GetPointDepth(float *pos);
/**
 * @brief GetPixelDepth get the depth of a pixel on screen
 * @param pixel_position_y top pixel has position 0
 */
float GetPixelDepth(int pixel_position_x, int pixel_position_y);
/**
 * @brief GetSelectionRay given a pixel, return a ray from camera to the pixel
 */
void GetSelectionRay(int pixel_position_x, int pixel_position_y, float *starting_point, float *ending_point);
void GetSelectionRay(int pixel_position_x, int pixel_position_y, double *starting_point, double *ending_point);
/**
 * @brief GetPixelWorldPosition given a pixel with depth value, return its world position
 */
void GetPixelWorldPosition(int pixel_position_x, int pixel_position_y, double depth, double *world_pos);
void GetPixelWorldPosition(int pixel_position_x, int pixel_position_y, float depth, float *world_pos);
void DisplayString(const char *str);
void IntersectWithHorizontalPlane(int pixel_position_x, int pixel_position_y, float &pos_x, float &pos_y);
void DrawSphere(float radius, int lats = 10, int longs = 10);
void DrawCylinder(float radius, float height, int slice = 10, int stack = 10);
void DrawCheckBoard(float width, float height, int x_slice, int y_slice, const float *color1, const float *color2, bool enable_lighting = false);
template <class FloatType>
inline void DrawArrow(FloatType start_x, FloatType start_y, FloatType start_z,
                      FloatType end_x, FloatType end_y, FloatType end_z,
                      FloatType arrow_head_length = 0.06,
                      FloatType arrow_head_radius = 0.02)
{
  const double kRadianPerDegree = 0.0174533;
  FloatType x = end_x - start_x;
  FloatType y = end_y - start_y;
  FloatType z = end_z - start_z;
  FloatType L = sqrt(x * x + y * y + z * z);
  glBegin(GL_LINES);
  glVertex3d(start_x, start_y, start_z);
  glVertex3d(end_x, end_y, end_z);
  glEnd();
  static GLUquadric *cyl = gluNewQuadric();
  glPushMatrix ();
  glTranslated(start_x, start_y, start_z);
  if ((x != 0.) || (y != 0.))
  {
    glRotated(atan2(y, x) / kRadianPerDegree, 0., 0., 1.);
    glRotated(atan2(sqrt(x * x + y * y), z) / kRadianPerDegree, 0., 1., 0.);
  }
  else if (z < 0)
  {
    glRotated(180, 1., 0., 0.);
  }
  glTranslatef(0, 0, L - arrow_head_length);
  gluQuadricDrawStyle (cyl, GLU_FILL);
  gluQuadricNormals (cyl, GLU_SMOOTH);
  gluCylinder(cyl, arrow_head_radius, 0.0, arrow_head_length, 12, 1);
  //  gluDeleteQuadric(cyl);
  glPopMatrix();

}

template <class FloatType>
inline void DrawArrow(FloatType *start, FloatType *end,
                      FloatType arrow_head_length = 0.06,
                      FloatType arrow_head_radius = 0.02)
{
  DrawArrow<FloatType>(start[0], start[1], start[2],
                       end[0], end[1], end[2],
                       arrow_head_length, arrow_head_radius);
}

template <class FloatType>
inline void DrawArrow(FloatType *start, FloatType *direction,  bool normlize_direction,
                      FloatType scale = 1.0,
                      FloatType arrow_head_length = 0.06,
                      FloatType arrow_head_radius = 0.02)
{
  FloatType end[3] = {direction[0], direction[1], direction[2]};
  if (normlize_direction)
  {
    float length = sqrtf(end[0] * end[0] + end[1] * end[1] + end[2] * end[2]);
    if (length > 1e-8)
    {
      float inv_length = FloatType(1.0) / length;
      end[0] *= inv_length;
      end[1] *= inv_length;
      end[2] *= inv_length;
    }
  }
  end[0] = start[0] + scale * end[0];
  end[1] = start[1] + scale * end[1];
  end[2] = start[2] + scale * end[2];
  DrawArrow<FloatType>(start[0], start[1], start[2],
                       end[0], end[1], end[2],
                       arrow_head_length, arrow_head_radius);
}

void DrawGradientBackGround(float *lower_color = NULL, float *upper_color = NULL);
void DrawAxis();


inline void Vertex3(float a, float b, float c) { glVertex3f(a, b, c); }
inline void Vertex3v(float* vertex) { glVertex3fv(vertex); }
inline void Vertex2(float* vertex) { glVertex2fv(vertex); }
inline void Normal(float* normal) { glNormal3fv(normal); }

inline void Vertex3(double a, double b, double c) { glVertex3d(a, b, c); }
inline void Vertex3v(double* vertex) { glVertex3dv(vertex); }
inline void Vertex2(double* vertex) { glVertex2dv(vertex); }
inline void Normal(double* normal) { glNormal3dv(normal); }


#endif // OPENGL_HELPER_H_
