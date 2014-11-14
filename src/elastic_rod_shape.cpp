#include "elastic_rod_shape.h"
#include "rainbow_color.h"
#include "opengl_helper.h"

ElasticRodShape::ElasticRodShape(int vert_num, real length) {
  // Vertex data
  vert_.resize(vert_num);
  // Edge data
  edge_vert_.resize(vert_num - 1);
  darboux_vector_.resize(vert_num - 2);

  real segment_length = length / (vert_num - 1);
  avg_edge_length_ = segment_length;
  Vec3 start(0.0, 1.0, 0.0);
//  Vec3 direction(segment_length, 0, 0);
  Vec3 direction(0, -segment_length, 0);
  Vec3 edge_direction(0, 0, segment_length);
  for (unsigned v = 0; v < vert_.size(); ++v) {
    vert_[v] = start + v * direction;
  }
  for (unsigned e = 0; e < edge_vert_.size(); ++e) {
    Vec3 mid = real(0.5) * (vert_[e] + vert_[e + 1]);
    edge_vert_[e] = mid + edge_direction;
  }
}

ElasticRodShape::ElasticRodShape(const ElasticRodShape &state) {
  (*this) = state;
}

ElasticRodShape &ElasticRodShape::operator=(const ElasticRodShape &state) {
  this->vert_ = state.vert_;
  this->edge_vert_ = state.edge_vert_;
  this->avg_edge_length_ = state.avg_edge_length_;
  this->darboux_vector_ = state.darboux_vector_;
  return *this;
}

void ElasticRodShape::Render() {
  // Render vertex
  glPointSize(4.0);
  glBegin(GL_POINTS);
  glColor3fv(kRed());
  for (auto v : vert_) {
    Vertex3(v[0], v[1], v[2]);
  }

  glColor3fv(kYellow());
  for (auto e : edge_vert_) {
    Vertex3(e[0], e[1], e[2]);
  }
  glEnd();

  // Render edge
  glColor3fv(kGreen());
  glBegin(GL_LINES);
  for (unsigned i = 0; i < edge_vert_.size(); ++i) {
    Vertex3(vert_[i][0], vert_[i][1], vert_[i][2]);
    Vertex3(vert_[i + 1][0], vert_[i + 1][1], vert_[i + 1][2]);
  }
  glEnd();
}
