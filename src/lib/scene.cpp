#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include "scene.h"
#include "config_file.h"
#include "macro_constant.h"
//#include "print_macro.h"

#define READ_VALUE(token) token = conf_->Get<decltype(token)>(STRINGIZE_TOKEN(token))
#define WRITE_VALUE(token) conf_->Put<decltype(token)>(STRINGIZE_TOKEN(token), token)

Scene::Scene(const char *scene_file_name)
//  : conf_(NULL)
{
  conf_ = NULL;
  Load(scene_file_name);
}

Scene::~Scene()
{
  delete conf_;
}
#include "vector_lib.h"
void Scene::Load(const char *scene_file_name)
{
  scene_file_name_ = std::string(scene_file_name);
  if (!conf_) delete conf_;
  conf_ = new ConfigFile(scene_file_name_.c_str());

  READ_VALUE(swing_angle);
  READ_VALUE(elevate_angle);
  READ_VALUE(screen_width);
  READ_VALUE(screen_height);
  READ_VALUE(zoom);
  READ_VALUE(field_of_view);
  READ_VALUE(near_plane);
  READ_VALUE(far_plane);
  READ_VALUE(zoom_step);
  READ_VALUE(min_zoom);
  READ_VALUE(target);

//  P(screen_width);
//  P(screen_height);
//  P(swing_angle);
//  P(elevate_angle);
//  P(zoom);
//  P(field_of_view);
//  P(near_plane);
//  P(far_plane);
//  P(zoom_step);
//  P(min_zoom);
//  P(dj::Vec3f(target));
}



void Scene::Save()
{
  WRITE_VALUE(screen_width);
  WRITE_VALUE(screen_width);
  WRITE_VALUE(swing_angle);
  WRITE_VALUE(elevate_angle);
  WRITE_VALUE(zoom);
  WRITE_VALUE(field_of_view);
  WRITE_VALUE(near_plane);
  WRITE_VALUE(far_plane);
  WRITE_VALUE(zoom_step);
  WRITE_VALUE(min_zoom);
  conf_->Save(scene_file_name_.c_str());
}

#undef READ_VALUE
#undef WRITE_VALUE
