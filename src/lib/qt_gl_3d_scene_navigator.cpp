#include "qt_gl_3d_scene_navigator.h"
#include <QPoint>
#include "scene.h"
#include "open_gl_qt.h"
#include "print_macro.h"
#include "global.h"

enum MouseState {
  kDefault,
  kZoom,
  kRotate	,
  kTranslate
};

int QtGL3DSceneNavigator::HandleMouseMove(QMouseEvent *event)
{
  QPoint pos = event->pos();
  switch (mouse_motion_mode_) {
    case kRotate:
      scene_->swing_angle   += (float)(pos.x() - last_mouse_pos_.x()) * 360 / (float)scene_->screen_width;
      scene_->elevate_angle += (float)(pos.y() - last_mouse_pos_.y()) * 180 / (float)scene_->screen_height;
      if  (scene_->elevate_angle > 90) {
        scene_->elevate_angle =  90;
      } else if (scene_->elevate_angle < -90) {
        scene_->elevate_angle = -90;
      }
      break;
    case kTranslate:
      scene_->target[0] -= (last_mouse_pos_.x() - pos.x()) * 0.02f;
      scene_->target[1] += (last_mouse_pos_.y() - pos.y()) * 0.02f;
      break;
    case kZoom:
      scene_->zoom -= (last_mouse_pos_.y() - pos.y()) * scene_->zoom_step;
      break;
    default:
      return -1;
      break;
  }
  last_mouse_pos_ = pos;
  return kPriority_;
  return -1;
}

int QtGL3DSceneNavigator::HandleMousePress(QMouseEvent* event)
{
  QPoint pos = event->pos();
  if (event->modifiers() & Qt::ShiftModifier) {
    mouse_motion_mode_ = kZoom;
  }
  //    else if (event->modifiers() & Qt::ControlModifier) {
  //      mouse_motion_mode_ = kTranslate;
  //    }
  else {
    mouse_motion_mode_ = kRotate;
  }
end:
  last_mouse_pos_ = pos;
  return kPriority_;
}

int QtGL3DSceneNavigator::HandleMouseRelease(QMouseEvent *event)
{
  QPoint pos = event->pos();
  mouse_motion_mode_ = kDefault;
  last_mouse_pos_ = pos;
  return -1;
}

int QtGL3DSceneNavigator::HandleKeyPress(QKeyEvent *event)
{
  const float kTranslateStep = 0.01f;
  switch (event->key()) {
    case Qt::Key_Q:
      exit(0);
      break;
    case Qt::Key_P:
      gl_widget_->ScreenShot("test.bmp", "bmp");
      break;
    case Qt::Key_W:
      scene_->Save();
      emit gl_widget_->updateStatusMessage("Scene file saved");
      break;
    case Qt::Key_A:
      scene_->zoom -= scene_->zoom_step;
      if (scene_->zoom < scene_->min_zoom) scene_->zoom = scene_->min_zoom;
      break;
    case Qt::Key_Z:
      scene_->zoom += scene_->zoom_step;
      break;
    case Qt::Key_X:
      draw_axis_ = !draw_axis_;
      break;
    case Qt::Key_Down:
      scene_->target[1] -= kTranslateStep;
      break;
    case Qt::Key_Up:
      scene_->target[1] += kTranslateStep;
      break;
    case Qt::Key_Left:
      scene_->target[0] -= kTranslateStep;
      break;
    case Qt::Key_Right:
      scene_->target[0] += kTranslateStep;
      break;
    case Qt::Key_Space:
      if (event->modifiers() & Qt::ShiftModifier) {
        global::pause_per_frame = !global::pause_per_frame;
      } else {
        global::simulate = !global::simulate;
      }
      break;
    default:
      break;
  }
  return -1;
}

int QtGL3DSceneNavigator::HandleKeyRelease(QKeyEvent *event)
{
  Q_UNUSED(event);
  return -1;
}

int QtGL3DSceneNavigator::PaintGL()
{
  if (draw_axis_) {
    DrawAxis();
  }
  return -1;
}

void QtGL3DSceneNavigator::set_gl_widget(OpenGLQt *widget)
{
  gl_widget_ = widget;
  scene_ = widget->scene();
}

void QtGL3DSceneNavigator::set_scene(Scene *scene) { scene_ = scene; }
