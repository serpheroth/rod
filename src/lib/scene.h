#ifndef SCENE_H
#define SCENE_H
#include <string>
class ConfigFile;
class Scene {
public:
  int		screen_width, screen_height;
  float	swing_angle, elevate_angle, zoom;
  float field_of_view, near_plane, far_plane;
  float zoom_step;
  float min_zoom;
  float* target;
  Scene() {}
  Scene(const char* scene_file_name);
  ~Scene();
  void Save(void);
  void Load(const char* scene_file_name);

//  template <class T>
//  T ReadValue(const char* key);

  template <class T>
  void ReadValue(const char* key, T*& value);

private:
  std::string scene_file_name_;
  ConfigFile* conf_;
}; // class Scene

#endif // SCENE_H
