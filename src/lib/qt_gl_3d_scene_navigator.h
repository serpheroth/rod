#ifndef QT_GL_3D_INPUT_HANLDER_H
#define QT_GL_3D_INPUT_HANLDER_H
#include "input_handler.h"
#include "singleton.h"

class OpenGLQt;
class Scene;
class QtGL3DSceneNavigator : public InputHandler {
public:
  // Mouse handlers
  virtual int HandleMouseMove(QMouseEvent* event);
  virtual int HandleMousePress(QMouseEvent *event);
  virtual int HandleMouseRelease(QMouseEvent* event);

  virtual int HandleKeyPress(QKeyEvent* event);
  virtual int HandleKeyRelease(QKeyEvent* event);

  virtual int PaintGL();
  void set_gl_widget(OpenGLQt* widget);
  void set_scene(Scene* scene);
private:
  bool draw_axis_;
  OpenGLQt* gl_widget_;
  Scene* scene_;
  QPoint last_mouse_pos_;
  QtGL3DSceneNavigator() : InputHandler(0), draw_axis_(true) {}
  int mouse_motion_mode_;// = kDefault;
  DECLARE_SINGLETON_CLASS(QtGL3DSceneNavigator)
};

#endif // QT_GL_3D_INPUT_HANLDER_H
