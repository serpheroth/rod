#include <iostream>
#include <fstream>
#include <cmath>
#include "opengl_helper.h"
//#include "opengl.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include "vector_lib.h"
#include "print_macro.h"

#if 0
void DisplayString(const char* str)
{
  glColor3f(0, 1, 0) ;
  glMatrixMode(GL_PROJECTION) ;
  glPushMatrix();
  glLoadIdentity();
  int screen_width = glutGet(GLUT_WINDOW_WIDTH);
  int screen_height = glutGet(GLUT_WINDOW_HEIGHT);
  gluOrtho2D(0, screen_width, 0, screen_height) ;
  glMatrixMode(GL_MODELVIEW) ;
  glPushMatrix();
  glLoadIdentity() ;
  glRasterPos2i(5, 5) ;
  for (int i = 0; str[i] != '\0'; i++) {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]) ;
  }
  glMatrixMode(GL_PROJECTION) ;
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW) ;
  glPopMatrix();
}
#endif

int WritePPM(const char* output_file_name, unsigned char* image, int resolution_x, int resolution_y)
{
  using namespace std;
  ofstream image_file;
  image_file.open(output_file_name) ;
  if (!image_file.is_open()) {
    std::cerr << "failed to open file for writing image." << std::endl ;
    return -1 ;
  }
  image_file << "P3" << std::endl ;
  image_file << resolution_x << " " << resolution_y << std::endl ;
  image_file << 255 << std::endl ;
  for (int y = 0 ; y < resolution_y ; y++) {
    for (int x = 0 ; x < resolution_x ; x++ ) {
      int pos = (y * resolution_x + x) * 3 ;
      image_file << ( (unsigned int) image[pos]) << " " << ( (unsigned int) image[pos + 1]) << " " << ((unsigned int) image[pos + 2]) << " " ;
    }
    image_file << std::endl ;
  }
  image_file.close() ;
  return  1 ;
}

double GetPointDepth(double* pos)
{
  double winX, winY, winZ; // holder for world coordinates
  GLint view_port[4]; // viewport dimensions+pos
  GLdouble projectioin_matrix[16];
  GLdouble model_view_matrix[16];
  glGetDoublev (GL_MODELVIEW_MATRIX, model_view_matrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projectioin_matrix);
  glGetIntegerv(GL_VIEWPORT, view_port);
  gluProject(pos[0], pos[1], pos[2],
             model_view_matrix, projectioin_matrix, view_port,
             &winX, &winY, &winZ);
  return winZ;
}

float GetPointDepth(float* pos)
{
  double tmp_pos[3] = {pos[0], pos[1], pos[2]};
  return (float) GetPointDepth(tmp_pos);
}

float GetPixelDepth(int pixel_position_x, int pixel_position_y)
{
  GLint view_port[4];
  glGetIntegerv(GL_VIEWPORT, view_port);
  float depth;
  glReadPixels(pixel_position_x, view_port[3] - pixel_position_y, 1, 1,
               GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
  return depth;
}


void GetPixelWorldPosition(int pixel_position_x, int pixel_position_y, double depth, double* world_pos)
{
  GLint view_port[4]; // viewport dimensions+pos
  GLdouble projectioin_matrix[16];
  GLdouble model_view_matrix[16];
  glGetDoublev (GL_MODELVIEW_MATRIX, model_view_matrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projectioin_matrix);
  glGetIntegerv(GL_VIEWPORT, view_port);
  //view[3]-cursorY = conversion from upper left (0,0) to lower left (0,0)
  //Unproject 2D Screen coordinates into wonderful world coordinates
  gluUnProject((GLdouble) pixel_position_x, (GLdouble) view_port[3] - pixel_position_y, depth,
               model_view_matrix, projectioin_matrix, view_port,
               world_pos + 0, world_pos + 1, world_pos + 2);
}

void GetPixelWorldPosition(int pixel_position_x, int pixel_position_y, float depth, float* world_pos) {
  double tmp_pos[3];
  GetPixelWorldPosition(pixel_position_x, pixel_position_y, (double) depth, tmp_pos);
  world_pos[0] = float(tmp_pos[0]);
  world_pos[1] = float(tmp_pos[1]);
  world_pos[2] = float(tmp_pos[2]);
}

void GetSelectionRay(int pixel_position_x, int pixel_position_y, double *starting_point, double *ending_point)
{
  double objX, objY, objZ; // holder for world coordinates
  GLint view_port[4]; // viewport dimensions+pos
  GLdouble projectioin_matrix[16];
  GLdouble model_view_matrix[16];
  glGetDoublev (GL_MODELVIEW_MATRIX, model_view_matrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projectioin_matrix);
  glGetIntegerv(GL_VIEWPORT, view_port);
  //view[3]-cursorY = conversion from upper left (0,0) to lower left (0,0)
  //Unproject 2D Screen coordinates into wonderful world coordinates
  gluUnProject((GLdouble) pixel_position_x, (GLdouble) view_port[3] - pixel_position_y, 1,
               model_view_matrix, projectioin_matrix, view_port, &objX, &objY, &objZ);
  ending_point[0] = objX;
  ending_point[1] = objY;
  ending_point[2] = objZ;
  gluUnProject((GLdouble) pixel_position_x, (GLdouble) view_port[3] - pixel_position_y, 0,
               model_view_matrix, projectioin_matrix, view_port, &objX, &objY, &objZ);
  starting_point[0] = objX;
  starting_point[1] = objY;
  starting_point[2] = objZ;
}

void GetSelectionRay(int pixel_position_x, int pixel_position_y, float *starting_point, float *ending_point)
{
  double start[3], end[3];
  GetSelectionRay(pixel_position_x, pixel_position_y, start, end);
  starting_point[0] = (float) start[0];
  starting_point[1] = (float) start[1];
  starting_point[2] = (float) start[2];
  ending_point[0] = (float) end[0];
  ending_point[1] = (float) end[1];
  ending_point[2] = (float) end[2];
}

void IntersectWithHorizontalPlane(int pixel_position_x, int pixel_position_y, float &pos_x, float &pos_y)
{
  float starting_point[3];
  float ending_point[3];
  GetSelectionRay(pixel_position_x, pixel_position_y, starting_point, ending_point);
  pos_x = -starting_point[1] / (ending_point[1] - starting_point[1]) * (ending_point[0] - starting_point[0]) + starting_point[0];
  pos_y = -starting_point[1] / (ending_point[1] - starting_point[1]) * (ending_point[2] - starting_point[2]) + starting_point[2];
}

void DrawSphere(float radius, int lats, int longs)
{

  GLUquadricObj* Sphere = gluNewQuadric();
  gluSphere(Sphere, radius, lats, longs);
  //  gluDeleteQuadric(Sphere);
  //  const float kPi = 3.1415926f;
  //  for (int i = 0; i <= lats; i++) {
  //    float lat0 = kPi * (-0.5 + (float) (i - 1) / lats);
  //    float z0  = radius * sin(lat0);
  //    float zr0 =  radius* cos(lat0);
  //
  //    float lat1 = kPi * (-0.5 + (float) i / lats);
  //    float z1 = radius * sin(lat1);
  //    float zr1 = radius * cos(lat1);
  //
  //    glBegin(GL_QUAD_STRIP);
  //    for (int j = 0; j <= longs; j++) {
  //      float lng = 2 * kPi * (float) (j - 1) / longs;
  //      float x = cos(lng);
  //      float y = sin(lng);
  //      glNormal3f(x * zr0, y * zr0, z0);
  //      glVertex3f(x * zr0, y * zr0, z0);
  //      glNormal3f(x * zr1, y * zr1, z1);
  //      glVertex3f(x * zr1, y * zr1, z1);
  //    }
  //    glEnd();
  //  }
}


void DrawCylinder(float radius, float height, int slice, int stack)
{
  static GLUquadricObj *quadObj = gluNewQuadric();
  gluCylinder(quadObj,
              radius, // base radius
              radius, // top radius
              height, // height
              slice, // slice
              stack // stack
             );
}


void DrawGradientBackGround(float *lower_color, float *upper_color) {
  // May background
//  float default_lower_color[3] = {30 / 255.0f, 30 / 255.0f, 30 / 255.0f};
//  float default_upper_color[3] = {130 / 255.0f, 140 / 255.0f, 160 / 255.0f};
  // Meshlab backgroud
  float default_lower_color[3] = {115 / 255.0f, 115 / 255.0f, 230 / 255.0f};
  float default_upper_color[3] = {0 / 255.0f, 0 / 255.0f, 0 / 255.0f};
  // White back ground
//  float default_lower_color[3] = {250 / 255.0f, 250 / 255.0f, 250 / 255.0f};
//  float default_upper_color[3] = {220 / 255.0f, 220 / 255.0f, 220 / 255.0f};

  if (lower_color == NULL) { lower_color = default_lower_color; }
  if (upper_color == NULL) { upper_color = default_upper_color; }
  GLboolean lighting_enabled, depth_test_enabled;
  glGetBooleanv(GL_DEPTH_TEST, &depth_test_enabled);
  glGetBooleanv(GL_LIGHTING, &lighting_enabled);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  glBegin(GL_QUADS);
  glDisable(GL_DEPTH_TEST);
  glColor3fv(lower_color);
  glVertex3f(-1.0,-1.0, 1);
  glVertex3f(1.0,-1.0, 1);
  glColor3fv(upper_color);
  glVertex3f(1.0, 1.0, 1);
  glVertex3f(-1.0, 1.0, 1);
  glEnd();
  if (lighting_enabled) { glEnable(GL_LIGHTING);}
  if (depth_test_enabled) { glEnable(GL_DEPTH_TEST);}

  glPopMatrix(); // Pop model view matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
}


void DrawAxis() {
  GLboolean lighting_enabled;
  glGetBooleanv(GL_LIGHTING, &lighting_enabled);
  glDisable(GL_LIGHTING);
  // X axis
  glColor3f(1, 0, 0) ;
  DrawArrow<float>(0.0, 0, 0, 1, 0, 0);
  // Y axis
  glColor3f(0, 1, 0) ;
  DrawArrow<float>(0.0, 0, 0, 0, 1, 0);
  // Z axis
  glColor3f(0, 0, 1);
  DrawArrow<float>(0.0, 0, 0, 0, 0, 1);
  glPopMatrix();
  if (lighting_enabled) { glEnable(GL_LIGHTING);}
}

///////////////////////////////////////////////////////////////////////////////
// write 2d text using GLUT
// The projection matrix must be set to orthogonal before call this function.
///////////////////////////////////////////////////////////////////////////////
void drawString(const char *str, int x, int y, float color[4], void *font)
{
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT); // lighting and color mask
    glDisable(GL_LIGHTING);     // need to disable lighting for proper text color
    glDisable(GL_TEXTURE_2D);

    glColor4fv(color);          // set text color
    glRasterPos2i(x, y);        // place text position

    // loop all characters in the string
    while(*str)
    {
        glutBitmapCharacter(font, *str);
        ++str;
    }

    glEnable(GL_TEXTURE_2D);
    glEnable(GL_LIGHTING);
    glPopAttrib();
}



///////////////////////////////////////////////////////////////////////////////
// draw a string in 3D space
///////////////////////////////////////////////////////////////////////////////
void drawString3D(const char *str, float pos[3], float color[4], void *font)
{
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT); // lighting and color mask
    glDisable(GL_LIGHTING);     // need to disable lighting for proper text color
    glDisable(GL_TEXTURE_2D);

    glColor4fv(color);          // set text color
    glRasterPos3fv(pos);        // place text position

    // loop all characters in the string
    while(*str)
    {
        glutBitmapCharacter(font, *str);
        ++str;
    }

    glEnable(GL_TEXTURE_2D);
    glEnable(GL_LIGHTING);
    glPopAttrib();
}

void DrawCheckBoard(float width, float height,
                    int x_slice, int y_slice,
                    const float *color1, const float *color2,
                    bool enable_lighting)
{
  float x_slice_size = width / x_slice;
  float y_slice_size = height / y_slice;
  glPushAttrib(GL_ENABLE_BIT);
  if (enable_lighting) {
    glEnable(GL_LIGHTING);
  }
  glEnable(GL_DEPTH_TEST);
  glBegin(GL_QUADS);
  float x = -width / 2;// + i * x_slice_size;
  for (int i = 0; i < x_slice + 1; ++i, x += x_slice_size) {
    float y = -height / 2;// + j * y_slice_size;
    for (int j = 0; j < y_slice + 1; ++j, y += y_slice_size) {
      const float* color = ((i + j) % 2 == 0) ? color1 : color2;
      if (enable_lighting) {
        float diffuse[4] = {color[0], color[1], color[2], 1};
//        float diffuse[4] = {0.4, 0.8, 0.2, 1};
        float ambient[4] = {0.2f, 0.2f, 0.2f, 1};
//        float specular[4] = {1, 1, 1, 1};
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
//        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
        glNormal3f(0, 0, 1);
      } else {
        glColor3fv(color);
      }
      glVertex2f(x, y);
      glVertex2f(x + x_slice_size, y);
      glVertex2f(x + x_slice_size, y + y_slice_size);
      glVertex2f(x, y + y_slice_size);
    }
  }
  glEnd();
  glPopAttrib();
}
