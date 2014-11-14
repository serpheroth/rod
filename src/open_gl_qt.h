#ifndef OPENGL_QT_H
#define OPENGL_QT_H
#include "opengl_helper.h"
#include <QGLWidget>
#include <vector>
#include <memory>

class Scene;
class InputHandler;

class OpenGLQt : public QGLWidget {
  Q_OBJECT
public:
  explicit OpenGLQt(QWidget *parent = 0);
  ~OpenGLQt();

  Scene* scene() { return scene_.get(); }
  void ScreenShot(const char* file_name, const char* image_type);
  QMouseEvent adjustPosForRetinaDisplay(QMouseEvent* e);

private:
  void keyReleaseEvent(QKeyEvent *event);
  void mousePressEvent(QMouseEvent *event);
  void mouseReleaseEvent(QMouseEvent *e);
  void mouseMoveEvent(QMouseEvent *event);
  void wheelEvent(QWheelEvent *e);
  void keyPressEvent(QKeyEvent *event);

signals:
  void updateStatusMessage(QString);
  void wheelUp();
  void wheelDown();

public slots:
  void idle();
  void saveScene();

protected:
  void initializeGL();
  void resizeGL(int w, int h);
  void paintGL();

private:
  std::auto_ptr<Scene> scene_;
  std::vector<InputHandler*> input_handlers_;
  QTimer* timer_;
};

#endif
