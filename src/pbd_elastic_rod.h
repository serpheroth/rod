#ifndef PBD_ELASTIC_ROD_H
#define PBD_ELASTIC_ROD_H
#include <vector>
#include "elastic_rod_shape.h"
#include "global.h"
#include "vector_lib.h"

class FrameAttachment;
class PBDElasticRod : public ElasticRodShape
{
public:
  typedef ElasticRodShape super;
  typedef dj::Vec<Vec3, 3, false> Mat3;
  PBDElasticRod(int vert_num, real length);
  ~PBDElasticRod();
  void Simulate(real dt);

  inline void EnforceEdgeConstraint(int e);
  inline void EnforceBisectorConstraint(int e);
  inline void EnforceGhostEdgeConstraint(int e);
  inline void EnforceBendingTwistingConstraint(int n);
  inline void ComputeOneDarbouxVector(ElasticRodShape *rod, int n, Vec3 &darboux_vector);
  void ComputeDarbouxVector(ElasticRodShape* rod, std::vector<Vec3>& darboux_vector);
  inline void ComputeMaterialFrameDerivative(int e,
                                             Mat3& d1p0, Mat3& d1p1, Mat3& d1p2,
                                             Mat3& d2p0, Mat3& d2p1, Mat3& d2p2,
                                             Mat3& d3p0, Mat3& d3p1, Mat3& d3p2);

  inline static void ComputeOneMaterialFrame(Vec3& p0, Vec3& p1, Vec3& pe, Mat3& frame)
  {
    frame[2] = p1 - p0;
    frame[2].Normalize();
    Vec3 tmp = pe - p0;
    dj::Cross3(frame[2](), tmp(), frame[1]());
    frame[1].Normalize();
    dj::Cross3(frame[1](), frame[2](), frame[0]());
  }

  inline void ComputeDarbouxGradient(Vec3& darboux_vector, real length,
                                     Mat3& da, Mat3& db,
                                     Mat3 (*dajpi)[3], Mat3 (*dbjpi)[3],
                                     Mat3& omega_pa,
                                     Mat3& omega_pb,
                                     Mat3& omega_pc,
                                     Mat3& omega_pd,
                                     Mat3& omega_pe
                                    );
  inline void CorrectGravityForGhostPoints(real dt);
  inline void FixEdge(int e);
  void ComputeMaterialFrame(ElasticRodShape* rod);
  void Render();

  ElasticRodShape rest_shape_;
  std::vector<Mat3> d_; // Material fame
  std::vector<real> mass_;
  Vec3 gravity_;
  std::vector<Vec3> prev_vm_;
  std::vector<Vec3> vm_;
  std::vector<Vec3> vel_vert_;

  std::vector<real> edge_length_;
  std::vector<Vec3> vel_edge_;
  std::vector<real> mid_edge_length_;
  FrameAttachment* boundary_[2];
  int v_num_;
  int e_num_;
  int iteration_;
};

// Return [v]
inline PBDElasticRod::Mat3 GetCrossProductMatrix(const ElasticRodShape::Vec3& v)
{
  PBDElasticRod::Mat3 mat;
  mat[0][0] = 0;     mat[0][1] = -v[2]; mat[0][2] = +v[1];
  mat[1][0] = +v[2]; mat[1][1] = 0;     mat[1][2] = -v[0];
  mat[2][0] = -v[1]; mat[2][1] = +v[0]; mat[2][2] = 0;
  return mat;
}

inline void GetCrossProductMatrix(const ElasticRodShape::Vec3& v, PBDElasticRod::Mat3& mat) {
  mat[0][0] = 0;     mat[0][1] = -v[2]; mat[0][2] = +v[1];
  mat[1][0] = +v[2]; mat[1][1] = 0;     mat[1][2] = -v[0];
  mat[2][0] = -v[1]; mat[2][1] = +v[0]; mat[2][2] = 0;
}
#endif // PBD_ELASTIC_ROD_H
