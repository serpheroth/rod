#ifndef ELASTIC_ROD_STATE_H
#define ELASTIC_ROD_STATE_H

#include <vector>
#include "global.h"
#include "vector_lib.h"

class ElasticRodShape {
public:
  typedef dj::Vec<real, 3, false> Vec3;
  ElasticRodShape(int vert_num, real length);
  ElasticRodShape(const ElasticRodShape& state);
  ElasticRodShape& operator=(const ElasticRodShape& state);

  void Render();
  real avg_edge_length_;
  std::vector<Vec3> vert_;
  std::vector<Vec3> edge_vert_;
  std::vector<Vec3> darboux_vector_;
};

#endif // ELASTIC_ROD_STATE_H
