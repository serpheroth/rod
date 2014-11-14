//#include <GL/glew.h>
//#include <QGLFunctions>
#include <QMenu>
#include <QTimer>
#include <QDebug>
#include <QFileDialog>
#include <QString>
#include <cmath>
#include <fstream>
#include <QByteArray>
#include <QIODevice>
#include <QFont>

#include <QLabel>
#include <QByteArray>
#include <QImageWriter>
#include <QString>
#include <QDebug>
#include <QMouseEvent>
#include <sstream>

#include "open_gl_qt.h"
#include "input_handler.h"
#include "qt_gl_3d_scene_navigator.h"
#include "scene.h"
#include "rainbow_color.h"
#include "fps_calculator.h"
#include "global.h"
#include "print_macro.h"
#include "fps_calculator.h"
#include "timer.h"
#include "simulator.h"

#ifndef CURRENT_DIRECTORY
#define CURRENT_DIRECTORY STRINGIZE_TOKEN(DATA_DIRECTORY)
#endif

FPSCalculator<Timer, float> fps_calculator;
const std::string kSceneFileName(std::string(global::kPWD) + "scene.conf");

unsigned int idle_count = 0;

enum MouseState {
  kDefault,
  kZoom,
  kRotate	,
  kTranslate,
};

OpenGLQt::OpenGLQt(QWidget *parent)
  : QGLWidget(parent)
  , scene_(new Scene(kSceneFileName.c_str()))
{
  //  input_handlers_.push_back(Singleton<SkeletonBuilder>::Instance());
  QtGL3DSceneNavigator* handler = Singleton<QtGL3DSceneNavigator>::Instance();
  handler->set_gl_widget(this);
  handler->set_scene(scene());
  input_handlers_.push_back(handler);
  input_handlers_.push_back(Singleton<Simulator>::Instance());
  timer_ = new QTimer(this);
  connect(timer_, SIGNAL(timeout()), this, SLOT(idle()));
  timer_->start(1);
}

OpenGLQt::~OpenGLQt()
{
}

void OpenGLQt::ScreenShot(const char *file_name, const char *image_type)
{
  QImage* image;
  QImageWriter* image_writer;
  image = new QImage(width(), height(), QImage::Format_RGB888);
  glReadPixels(0, 0, width(), height(), GL_RGB, GL_UNSIGNED_BYTE, image->bits());
  //  image_writer = new QImageWriter("test.bmp", "bmp");
  image_writer = new QImageWriter(file_name, image_type);
  image_writer->setQuality(100);
  image_writer->write(image->mirrored(false, true));
  P(file_name, image_type);
  delete image;
  delete image_writer;
}

void OpenGLQt::saveScene()
{
  scene_->Save();
  emit updateStatusMessage(QString("Scene file save"));
}

void OpenGLQt::initializeGL(void)
{
  //  int arg = 0;
  //  glutInit(&arg, NULL);
  //  makeCurrent();
  ASSERT(glewInit() == GLEW_OK);
  glShadeModel(GL_SMOOTH);
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  //  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glClearDepth(1.0f);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  glEnable(GL_NORMALIZE);
}

void OpenGLQt::resizeGL(int w, int h)
{
  scene_->screen_width = w, scene_->screen_height = h;
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // (Field of view, aspect ratio, near plance, far plane)
  gluPerspective(scene_->field_of_view,  float(w) / h, scene_->near_plane, scene_->far_plane);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glEnable(GL_DEPTH_TEST);
  updateGL();
}

void OpenGLQt::paintGL()
{
  glClearColor(1, 1, 1, 1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  DrawGradientBackGround();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
#if 1
  {
    float LightDiffuse[] = {1.0, 1.0, 1.0, 0};
    float LightPosition[] = {0, 0, 0.00100, 0};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, LightDiffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);
    glEnable(GL_LIGHT0);
  }
#endif
  // (eye_pos, focus_point, view_up)
  gluLookAt(0, 0, scene_->zoom, 0, 0, 0, 0, 1, 0);
  glRotatef(scene_->elevate_angle, 1, 0, 0) ;
  glRotatef(scene_->swing_angle, 0, 1, 0) ;
  glTranslatef(-scene_->target[0], -scene_->target[1], -scene_->target[2]);
  for (InputHandler * handler : input_handlers_) {
    handler->PaintGL();
  }
  //  return;
  //------------------------------------------------------------------------------
  // Render texts on screen
  QFont sansFont("Helvetica [Cronyx]", 14);
  sansFont.setFamily("sans serif");
  int font_height = QFontMetrics(sansFont).height();
  glColor3fv(kBlack());
  int base = height() - 5;
  QString texts[] = {
    QString("FPS: %1").arg(fps_calculator.fps()),
    QString("Simulating: ") + ((global::simulate) ?  QString("yes") : QString("no")),
    QString("Pause every frame: ") + ((global::pause_per_frame) ?  QString("yes") : QString("no")),
  };
  for (unsigned i = 0; i < sizeof(texts) / sizeof(QString); ++i) {
    renderText(5, base - font_height * i, texts[i], sansFont);
  }
}

void OpenGLQt::keyPressEvent(QKeyEvent *event)
{
  for (auto handler : input_handlers_) {
    handler->HandleKeyPress(event);
  }
  updateGL();
}

void OpenGLQt::mousePressEvent(QMouseEvent *e)
{
  QMouseEvent event = adjustPosForRetinaDisplay(e);
  for (InputHandler * handler : input_handlers_) {
    handler->HandleMousePress(&event);
  }
}

void OpenGLQt::mouseReleaseEvent(QMouseEvent *e)
{
  QMouseEvent event = adjustPosForRetinaDisplay(e);
  for (InputHandler * handler : input_handlers_) {
    handler->HandleMouseRelease(&event);
  }
}

void OpenGLQt::mouseMoveEvent(QMouseEvent *e)
{
  QMouseEvent event = adjustPosForRetinaDisplay(e);
  for (InputHandler * handler : input_handlers_) {
    handler->HandleMouseMove(&event);
  }
  QPoint pos = event.pos();
  emit updateStatusMessage(QString("Cursor position: [%1, %2]").arg(pos.x()).arg(pos.y()));
  updateGL();
}

void OpenGLQt::wheelEvent(QWheelEvent *e)
{
  if (e->delta() > 0) {
    // wheel goes up
    scene_->zoom -= scene_->zoom_step;
    if (scene_->zoom < scene_->min_zoom) scene_->zoom = scene_->min_zoom;
  } else {
    scene_->zoom += scene_->zoom_step;
  }
  updateGL();
}

QMouseEvent OpenGLQt::adjustPosForRetinaDisplay(QMouseEvent *e)
{
  QPoint pos = e->pos();
  GLint view_port[4]; // viewport dimensions+pos
  glGetIntegerv(GL_VIEWPORT, view_port);
  pos *= view_port[2] / width();
  QMouseEvent tmp_event(e->type(), pos, e->button(), e->buttons(), e->modifiers());
  return tmp_event;
}

void OpenGLQt::idle()
{
  if (global::simulate) {
    for (InputHandler * handler : input_handlers_) {
      handler->Idle();
    }
    if (global::pause_per_frame) {
      global::simulate = false;
    }
  }
  fps_calculator.CalculateFPS();
  updateGL();
}

void OpenGLQt::keyReleaseEvent(QKeyEvent *event)
{
  for (InputHandler * handler : input_handlers_) {
    handler->HandleKeyRelease(event);
  }
}
