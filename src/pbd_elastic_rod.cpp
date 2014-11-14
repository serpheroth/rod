#include "opengl_helper.h"
#include "pbd_elastic_rod.h"
#include "rainbow_color.h"
#include "frame_attachment.h"

inline PBDElasticRod::Mat3
GetRotationMatrix(ElasticRodShape::Vec3& vec, real angle_in_radian)
{
  PBDElasticRod::Mat3 m;
  real cosine = cos(angle_in_radian);
  real sine = sin(angle_in_radian);
  m[0][0] = cosine + vec[0] * vec[0] * (1 - cosine);
  m[0][1] = vec[0] * vec[1] * (1 - cosine) - vec[2] * sine;
  m[0][2] = vec[0] * vec[2] * (1 - cosine) + vec[1] * sine;

  m[1][0] = vec[0] * vec[1] * (1 - cosine) + vec[2] * sine;
  m[1][1] = cosine + vec[1] * vec[1] * (1 - cosine);
  m[1][2] = vec[1] * vec[2] * (1 - cosine) - vec[0] * sine;

  m[2][0] = vec[0] * vec[2] * (1 - cosine) - vec[1] * sine;
  m[2][1] = vec[1] * vec[2] * (1 - cosine) + vec[0] * sine;
  m[2][2] = cosine + vec[2] * vec[2] * (1 - cosine);
  return m;
}



PBDElasticRod::PBDElasticRod(int vert_num, real length)
  : ElasticRodShape(vert_num, length)
  , rest_shape_(vert_num, length)
  , gravity_(global::gravity)
{
  // Vertex data
  v_num_ = vert_num;
  vel_vert_.resize(v_num_);

  // Edge data
  e_num_ = vert_num - 1;
  d_.resize(e_num_);
  vel_edge_.resize(e_num_);
  prev_vm_ = vm_ = vel_edge_;
  edge_length_.resize(e_num_);
  mid_edge_length_.resize(v_num_ - 2);
  for (int e = 0; e < e_num_; ++e) {
    edge_length_[e] = (rest_shape_.vert_[e + 1] - rest_shape_.vert_[e]).Magnitude();
  }
  for (int e = 0; e < e_num_ - 1; ++e) {
    mid_edge_length_[e] = real(0.5) * (edge_length_[e] + edge_length_[e + 1]);
    //    mid_edge_length_[e] = avg_edge_length_;
  }
  ComputeDarbouxVector(&rest_shape_, rest_shape_.darboux_vector_);
  ComputeDarbouxVector(this, darboux_vector_);
  Mat3 frame;
  ComputeOneMaterialFrame(vert_[0], vert_[1], edge_vert_[0], frame);
  boundary_[0] = new FrameAttachment(this, frame, frame, vert_[0]);

  ComputeOneMaterialFrame(vert_.back(), vert_[vert_.size() - 2], edge_vert_.back(), frame);
  boundary_[1] = new FrameAttachment(this, frame, frame, vert_.back());
}

PBDElasticRod::~PBDElasticRod()
{
  delete boundary_[0];
  delete boundary_[1];
}

void PBDElasticRod::Simulate(real dt)
{
  // Apply gravity
  OMP_FOR
  for (unsigned i =  0; i < vel_vert_.size(); ++i) {
    vel_vert_[i] = dt * gravity_;
  }
  CorrectGravityForGhostPoints(dt);
  OMP_FOR
  for (unsigned i =  0; i < vel_edge_.size(); ++i) {
    vel_edge_[i]  += dt * gravity_;
  }

  OMP_FOR
  for (unsigned i =  0; i < vel_vert_.size(); ++i) {
    vert_[i] += dt * vel_vert_[i];
  }
  OMP_FOR
  for (unsigned i =  0; i < vel_edge_.size(); ++i) {
    edge_vert_[i] += dt * vel_edge_[i];
  }

  std::vector<Vec3> tmp_vert_ = vert_;
  std::vector<Vec3> tmp_edge_vert_ = edge_vert_;
  Mat3 d;
  ComputeOneMaterialFrame(tmp_vert_[v_num_ - 2],
                          tmp_vert_[v_num_ - 1],
                          tmp_edge_vert_.back(), d);
  //    const real angle_in_radian = dj::Degree2Radian(90.0);
  //  const real angle_in_radian = 0.00001;
  //  Mat3 rotation_matrix = GetRotationMatrix(d[2], angle_in_radian);
  //  Vec3 new_angle_direction(rotation_matrix[0] * d[0], rotation_matrix[1] * d[0], rotation_matrix[2] * d[0]);
  //  Vec3 tip_edge = real(0.5) * (tmp_vert_[v_num_ - 2] + tmp_vert_[v_num_ - 1]) + avg_edge_length_ * new_angle_direction;
  //  edge_vert_.back() = tip_edge;
  static int i = 1;
  if (i == 0) {
    int test = 4;
    switch (test) {
      case 0:
        vert_.back()[0] = avg_edge_length_;
        vert_.back()[2] = +avg_edge_length_;
        edge_vert_.back()[0] = 0;//0.5 * avg_edge_length_;
        edge_vert_.back()[2] = +0.5 * avg_edge_length_;
        break;
      case 1:
        vert_.back()[0] = avg_edge_length_;
        vert_.back()[2] = -avg_edge_length_;
        edge_vert_.back()[0] = 2.0 * avg_edge_length_;
        edge_vert_.back()[2] = -0.5 * avg_edge_length_;
        break;
      case 2:
        vert_.back()[0] = avg_edge_length_;
        vert_.back()[1] += +avg_edge_length_;
        edge_vert_.back()[0] = 1.0 * avg_edge_length_;
        edge_vert_.back()[1] += +0.5 * avg_edge_length_;
        break;
      case 3:
        vert_.back()[0] = avg_edge_length_;
        vert_.back()[1] -= +avg_edge_length_;
        edge_vert_.back()[0] = 1.0 * avg_edge_length_;
        edge_vert_.back()[1] -= +0.5 * avg_edge_length_;
        break;
      case 4:
//        edge_vert_.back() = tip_edge;
        break;
      default:
        break;
    }
    ++i;
  }

  FrameAttachment::FrameUpdator frame_updator = [&](FrameAttachment & attachment) {
    static Vec3 axis(0, -1, 0);
    const real kAngularVelInDegree = 40;
    const real angle_in_radian = dj::Degree2Radian(kAngularVelInDegree * dt);
    Mat3 rot = GetRotationMatrix(axis, angle_in_radian);
    Mat3 tmp = attachment.current_frame_;
    dj::MulMatrixRightTransposed3x3<real>((real (*)[3]) &tmp[0][0], (real (*)[3]) &rot[0][0],
                                          (real (*)[3]) &attachment.current_frame_[0][0]);
  };
  boundary_[0]->UpdateFrame(frame_updator);

  FrameAttachment::FrameUpdator frame_updator1 = [&](FrameAttachment & attachment) {
    Vec3 axis = attachment.current_frame_[2];
    const real kAngularVelInDegree = 40;
    const real angle_in_radian = dj::Degree2Radian(kAngularVelInDegree * dt);
    Mat3 rot = GetRotationMatrix(axis, angle_in_radian);
    Mat3 tmp = attachment.current_frame_;
    dj::MulMatrixRightTransposed3x3<real>((real (*)[3]) &tmp[0][0], (real (*)[3]) &rot[0][0],
                                          (real (*)[3]) &attachment.current_frame_[0][0]);
  };
  boundary_[1]->UpdateFrame(frame_updator1);

  for (int i = 0; i < rod::pbd_iteration; ++i) {
    //#define SQUENTIAL_GAUSS_SEIDEL
#ifdef SQUENTIAL_GAUSS_SEIDEL
    for (int e = 0; e < e_num_; ++e) {
      EnforceEdgeConstraint(e);
      EnforceBisectorConstraint(e);
      EnforceGhostEdgeConstraint(e);
    }
    int n = 0;
    for (; n < v_num_ - 2; n += 1) {
      EnforceBendingTwistingConstraint(n);
    }
#else
    int left = 0, right = e_num_ - 1;
    for (; left < right; ++left, --right) {
      EnforceEdgeConstraint(left);
      EnforceBisectorConstraint(left);
      EnforceGhostEdgeConstraint(left);

      EnforceEdgeConstraint(right);
      EnforceBisectorConstraint(right);
      EnforceGhostEdgeConstraint(right);
    }
    if (left == right) {
      EnforceEdgeConstraint(left);
      EnforceBisectorConstraint(left);
      EnforceGhostEdgeConstraint(left);
    }

    left = 0, right = v_num_ - 3;
    for (; left < right; ++left, --right) {
      EnforceBendingTwistingConstraint(left);
      EnforceBendingTwistingConstraint(right);
    }
    if (left == right) {
      EnforceBendingTwistingConstraint(left);
    }
#endif
    boundary_[0]->EnforceFrameAttachment(vert_[0], vert_[1], edge_vert_[0], edge_length_[0]);
    boundary_[0]->EnforcePosAttachment(vert_[0]);
    boundary_[1]->EnforceFrameAttachment(vert_.back(), vert_[vert_.size() - 2], edge_vert_.back(), edge_length_.back());
    //    FixEdge(0);
    //    edge_vert_.back() = tip_edge;
    //    vert_[v_num_ - 2] = tmp_vert_[v_num_ - 2];
    //    vert_[v_num_ - 1] = tmp_vert_[v_num_ - 1];
  }

  for (unsigned i =  0; i < vel_vert_.size(); ++i) {
    vel_vert_[i] = (vert_[i] - tmp_vert_[i]) / dt;
  }
  for (unsigned i =  0; i < vel_edge_.size(); ++i) {
    vel_edge_[i] = (edge_vert_[i] - tmp_edge_vert_[i]) / dt;
  }
}

inline void PBDElasticRod::EnforceEdgeConstraint(int e)
{
  Vec3 direction = vert_[e + 1] - vert_[e];
  real length = direction.Magnitude();
  real scale = edge_length_[e] / length;
  scale *= real(0.5);
  direction *= scale;
  Vec3 center = real(0.5) * (vert_[e + 1] + vert_[e]);
  vert_[e + 1] = center + direction;
  vert_[e + 0] = center - direction;
}


inline void PBDElasticRod::EnforceBisectorConstraint(int e)
{
  Vec3& p0 = vert_[e];
  Vec3& p1 = vert_[e + 1];
  Vec3& p2 = edge_vert_[e];
  Vec3 p0p2 = p0 - p2;
  Vec3 p2p1 = p2 - p1;
  Vec3 p1p0 = p1 - p0;
  Vec3 pm = real(0.5f) * (p0 + p1);
  real lambda = p0p2.MagnitudeSquare() + p2p1.MagnitudeSquare() + p1p0.MagnitudeSquare();
  lambda = (p2 - pm) * p1p0 / lambda;
  p0 += -lambda * p0p2;
  p1 += -lambda * p2p1;
  p2 += -lambda * p1p0;
}

;

inline void PBDElasticRod::EnforceGhostEdgeConstraint(int e)
{
  Vec3& p0 = vert_[e];
  Vec3& p1 = vert_[e + 1];
  Vec3& p2 = edge_vert_[e];
  Vec3 pm = real(0.5f) * (p0 + p1);
  Vec3 p2pm = p2 - pm;
  real p2pm_mag = p2pm.Magnitude();
  p2pm *= (1 / p2pm_mag);
  real lambda = (p2pm_mag - avg_edge_length_) / real(0.25 + 0.25 + 1);
  p0 += +real(0.5) * lambda * p2pm;
  p1 += +real(0.5) * lambda * p2pm;
  p2 += -lambda * p2pm;
}

// n \in [0, v_num - 2)
inline void PBDElasticRod::EnforceBendingTwistingConstraint(int n)
{
  Vec3 darboux_vector;
  ComputeOneDarbouxVector(this, n, darboux_vector);

  {
    //    Vec3 grad;
    //    Vec3 tmp_darb = darboux_vector;
    //    int c = 2;
    //    //    Vec3& target = vert_.back();
    //    Vec3& target = edge_vert_.back();
    //    Vec3 cur = target;
    //    real diff = 1e-10;
    //    for (int i = 0; i < 3; ++i) {
    //      target = cur;
    //      target[i] += diff;
    //      ComputeOneDarbouxVector(this, n, darboux_vector);
    //      grad[i] = (darboux_vector[c] - tmp_darb[c]) / diff;
    //    }
    //    KK;
    //    P(grad);
    //    darboux_vector = tmp_darb;
    //    target = cur;
  }
  Vec3 constraint_value(rod::bending_twising_stiffness[0] * (darboux_vector[0] - rest_shape_.darboux_vector_[n][0]),
                        rod::bending_twising_stiffness[1] * (darboux_vector[1] - rest_shape_.darboux_vector_[n][1]),
                        rod::bending_twising_stiffness[2] * (darboux_vector[2] - rest_shape_.darboux_vector_[n][2]));
  {
    //  P(constraint_value);
    //  P(darboux_vector);
    //  ASSERT(darboux_vector == darboux_vector,
    //         P(n)
    //         P(vert_[n], vert_[n + 1], edge_vert_[n ])
    //         P(vert_[n + 1], vert_[n + 2], edge_vert_[n + 1])
    //        );
    //  if (constraint_value.Magnitude() < real(1e-6)) {
    //    return;
    //  }
    //  P(darboux_vector);
  }
  Mat3& da = d_[n];
  Mat3& db = d_[n + 1];
  ComputeOneMaterialFrame(vert_[n], vert_[n + 1], edge_vert_[n], da);
  ComputeOneMaterialFrame(vert_[n + 1], vert_[n + 2], edge_vert_[n + 1], db);
  {
    //      P(vert_[n], vert_[n + 1], edge_vert_[n]);
    //      P(vert_[n + 1], vert_[n + 2], edge_vert_[n + 1]);
    //      P(avg_edge_length_);
    //      P(da);
    //      P(db);
    //  P(da, db);
  }
  Mat3 dajpi[3][3];
  ComputeMaterialFrameDerivative(n,
                                 dajpi[0][0], dajpi[0][1], dajpi[0][2],
                                 dajpi[1][0], dajpi[1][1], dajpi[1][2],
                                 dajpi[2][0], dajpi[2][1], dajpi[2][2]);
  Mat3 dbjpi[3][3];
  ComputeMaterialFrameDerivative(n + 1,
                                 dbjpi[0][0], dbjpi[0][1], dbjpi[0][2],
                                 dbjpi[1][0], dbjpi[1][1], dbjpi[1][2],
                                 dbjpi[2][0], dbjpi[2][1], dbjpi[2][2]);
  {
    //    for (int i = 0; i < 3; ++i) {
    //      for (int j = 0; j < 3; ++j) {
    //        P(i, j); KK;
    //        P(dajpi[i][j]);
    //        P(dbjpi[i][j]);
    //      }
    //    }
    //        P(dbjpi[0][2]);
    //        P(dbjpi[1][2]);
    //        P(dbjpi[2][2]);
  }
  Mat3 constraint_jacobian[5];
  ComputeDarbouxGradient(darboux_vector, mid_edge_length_[n], da, db, dajpi, dbjpi,
                         constraint_jacobian[0],
                         constraint_jacobian[1],
                         constraint_jacobian[2],
                         constraint_jacobian[3],
                         constraint_jacobian[4]);
  Vec3* verts[] = {
    &vert_[n],
    &vert_[n + 1],
    &vert_[n + 2],
    &edge_vert_[n],
    &edge_vert_[n + 1]
  };
  //  P(constraint_value);
#if 1
#if 0
  for (int c = 0; c < 3; ++c) {
    for (int i = 0; i < 5; ++i) {
      *(verts[i]) -= (1e-8 * constraint_value[c] * constraint_jacobian[i][c]);
    }
  }
  //  P(constraint_value, darboux_vector);
#else
  for (int c = 0; c < 3; ++c) {
    real factor = 0;
    for (int i = 0; i < 5; ++i) {
      factor += constraint_jacobian[i][c] * constraint_jacobian[i][c];
    }
    //    P(1 / factor);
    factor = -constraint_value[c] / factor;
    factor *= 6e-2f;
    for (int i = 0; i < 5; ++i) {
      *(verts[i]) += (factor * constraint_jacobian[i][c]);
      //      ASSERT((*verts[i])[0] == (*verts[i])[0], P(factor, constraint_jacobian[i][c], *verts[i]));
      //      ASSERT((*verts[i])[1] == (*verts[i])[1], P(factor, constraint_jacobian[i][c], *verts[i]));
      //      ASSERT((*verts[i])[2] == (*verts[i])[2], P(factor, constraint_jacobian[i][c], *verts[i]));
    }
  }
#endif
#else
  Mat3 factor_matrix;
  Mat3 tmp_mat;
  for (int i = 0; i < 5; ++i) {
    dj::MulMatrixRightTransposed3x3<real>((real (*)[3]) &constraint_jacobian[i][0][0],
                                          (real (*)[3]) &constraint_jacobian[i][0][0],
                                          (real (*)[3]) &tmp_mat[0][0]);
    factor_matrix += tmp_mat;
  }
  Vec3 tmp(constraint_value[0], constraint_value[1], constraint_value[2]);
  dj::GaussianElimination<real, 3>(&factor_matrix[0][0], &tmp[0]);

  for (int i = 0; i < 5; ++i) {
    (*verts[i])[0] -= (constraint_jacobian[i][0][0] * tmp[0] + constraint_jacobian[i][1][0] * tmp[1] + constraint_jacobian[i][2][0] * tmp[2]);
    (*verts[i])[1] -= (constraint_jacobian[i][0][1] * tmp[0] + constraint_jacobian[i][1][1] * tmp[1] + constraint_jacobian[i][2][1] * tmp[2]);
    (*verts[i])[2] -= (constraint_jacobian[i][0][2] * tmp[0] + constraint_jacobian[i][1][2] * tmp[1] + constraint_jacobian[i][2][2] * tmp[2]);
    //    ASSERT(dj::Abs((*verts[i])[0]) < 4, P(*verts[i], constraint_jacobian[i], tmp));
    //    ASSERT(dj::Abs((*verts[i])[1]) < 4, P(*verts[i], constraint_jacobian[i], tmp));
    //    ASSERT(dj::Abs((*verts[i])[2]) < 4, P(*verts[i], constraint_jacobian[i], tmp));
  }
#endif
}


void PBDElasticRod::ComputeOneDarbouxVector(ElasticRodShape * rod, int n, Vec3 & darboux_vector)
{
  Mat3 d0, d1;
  ComputeOneMaterialFrame(rod->vert_[n], rod->vert_[n + 1], rod->edge_vert_[n], d0);
  ComputeOneMaterialFrame(rod->vert_[n + 1], rod->vert_[n + 2], rod->edge_vert_[n + 1], d1);
  //  P(n);
  //  P(rod->vert_[n], rod->vert_[n + 1], rod->edge_vert_[n]);
  //  P(rod->vert_[n + 1], rod->vert_[n + 2], rod->edge_vert_[n + 1]);
  //  P(d0);
  //  P(d1);
  real factor = 1 + d0[0] * d1[0] + d0[1] * d1[1] + d0[2] * d1[2];
  //  P(d0[0] * d1[0])P(d0[1] * d1[1])P(d0[2] * d1[2]);
  factor = real(2.0) / (mid_edge_length_[n] * factor);
  {
    //        int i = 0, j = 1, k = 2;
    int i = 0, j = 2, k = 1;
    darboux_vector[i] = d0[j] * d1[k] - d0[k] * d1[j];
  }
  {
    //        int i = 1, j = 2, k = 0;
    int i = 1, j = 0, k = 2;
    darboux_vector[i] = d0[j] * d1[k] - d0[k] * d1[j];
  }
  {
    //        int i = 2, j = 0, k = 1;
    int i = 2, j = 1, k = 0;
    darboux_vector[i] = d0[j] * d1[k] - d0[k] * d1[j];
  }
  darboux_vector *= factor;
}

void PBDElasticRod::ComputeDarbouxVector(ElasticRodShape * rod, std::vector<Vec3> &darboux_vector)
{
  ComputeMaterialFrame(rod);
  for (int e = 0; e < e_num_ - 1; ++e) {
    real factor = 1 + d_[e][0] * d_[e + 1][0] + d_[e][1] * d_[e + 1][1] + d_[e][2] * d_[e + 1][2];
    factor = real(2.0) / (mid_edge_length_[e] * factor);
    {
      //            int i = 0, j = 1, k = 2;
      int i = 0, j = 2, k = 1;
      darboux_vector[e][i] = d_[e][j] * d_[e + 1][k] - d_[e][k] * d_[e + 1][j];
    }
    {
      //            int i = 1, j = 2, k = 0;
      int i = 1, j = 0, k = 2;
      darboux_vector[e][i] = d_[e][j] * d_[e + 1][k] - d_[e][k] * d_[e + 1][j];
    }
    {
      //            int i = 2, j = 0, k = 1;
      int i = 2, j = 1, k = 0;
      darboux_vector[e][i] = d_[e][j] * d_[e + 1][k] - d_[e][k] * d_[e + 1][j];
    }
    darboux_vector[e] *= factor;
  }
}

void PBDElasticRod::ComputeMaterialFrameDerivative(int e,
                                                   Mat3 & d1p0, Mat3 & d1p1, Mat3 & d1p2,
                                                   Mat3 & d2p0, Mat3 & d2p1, Mat3 & d2p2,
                                                   Mat3 & d3p0, Mat3 & d3p1, Mat3 & d3p2)
{
  Vec3& p0 = vert_[e];
  Vec3& p1 = vert_[e + 1];
  Vec3& p2 = edge_vert_[e];
  //  P(p0);
  //  P(p1);
  //  P(p2);
  Mat3& d = d_[e];
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // d3pi
  Vec3 p01 = p1 - p0;
  real length_p01 = p01.Magnitude();
  d3p0[0] = d[2][0] * d[2];
  d3p0[1] = d[2][1] * d[2];
  d3p0[2] = d[2][2] * d[2];
  d3p0[0][0] -= 1;
  d3p0[1][1] -= 1;
  d3p0[2][2] -= 1;
  d3p0[0] *= (1 / length_p01);
  d3p0[1] *= (1 / length_p01);
  d3p0[2] *= (1 / length_p01);

  d3p1[0] = real(-1.0) * d3p0[0];
  d3p1[1] = real(-1.0) * d3p0[1];
  d3p1[2] = real(-1.0) * d3p0[2];

  d3p2[0].Fill(real(0.0));
  d3p2[0].Fill(real(0.0));
  d3p2[0].Fill(real(0.0));
  //  P(d3p0);
  //  P(d3p1);
  //  P(d3p2);
  //  exit(0);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // d2pi
  Vec3 p02 = p2 - p0;
  Vec3 p01_cross_p02;
  dj::Cross3(p01(), p02(), p01_cross_p02());
  real length_cross = p01_cross_p02.Magnitude();
  Mat3 mat(d[1][0] * d[1],
           d[1][1] * d[1],
           d[1][2] * d[1]);
  mat[0][0] -= 1;
  mat[1][1] -= 1;
  mat[2][2] -= 1;
  mat[0] *= (real(-1.0) / length_cross);
  mat[1] *= (real(-1.0) / length_cross);
  mat[2] *= (real(-1.0) / length_cross);
  //  Vec3 p2p1 = p2 - p1;
  Mat3 product_matrix;
  GetCrossProductMatrix(p2 - p1, product_matrix);
  dj::MulMatrix3x3<real>(&mat[0][0], &product_matrix[0][0], &d2p0[0][0]);
  GetCrossProductMatrix(p0 - p2, product_matrix);
  dj::MulMatrix3x3<real>(&mat[0][0], &product_matrix[0][0], &d2p1[0][0]);
  GetCrossProductMatrix(p1 - p0, product_matrix);
  //  P(p1);
  //  P(p0);
  //  P(p1 - p0);
  //  P(product_matrix);
  //  P(mat);
  //  exit(0);
  dj::MulMatrix3x3<real>(&mat[0][0], &product_matrix[0][0], &d2p2[0][0]);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // d1pi
  Mat3 product_mat_d3 = GetCrossProductMatrix(d[2]);
  Mat3 product_mat_d2 = GetCrossProductMatrix(d[1]);
  Mat3 m1, m2;

  dj::MulMatrix3x3<real>(&product_mat_d3[0][0], &d2p0[0][0], &m1[0][0]);
  dj::MulMatrix3x3<real>(&product_mat_d2[0][0], &d3p0[0][0], &m2[0][0]);
  d1p0 = m2 - m1;

  dj::MulMatrix3x3<real>(&product_mat_d3[0][0], &d2p1[0][0], &m1[0][0]);
  dj::MulMatrix3x3<real>(&product_mat_d2[0][0], &d3p1[0][0], &m2[0][0]);
  d1p1 = m2 - m1;

  dj::MulMatrix3x3<real>(&product_mat_d3[0][0], &d2p2[0][0], &d1p2[0][0]);
  d1p2[0] *= real(-1);
  d1p2[1] *= real(-1);
  d1p2[2] *= real(-1);
}



inline void PBDElasticRod::ComputeDarbouxGradient(Vec3& darboux_vector, real length,
                                                  Mat3& da, Mat3& db,
                                                  Mat3 (*dajpi)[3], Mat3 (*dbjpi)[3],
                                                  Mat3& omega_pa,
                                                  Mat3& omega_pb,
                                                  Mat3& omega_pc,
                                                  Mat3& omega_pd,
                                                  Mat3& omega_pe
                                                 )
{
  const int permutation[3][3] = {
    //0, 1, 2,
    //1, 2, 0,
    //2, 1, 0
    0, 2, 1,
    1, 0, 2,
    2, 1, 0
  };
  real x = 1 + da[0] * db[0] + da[1] * db[1] + da[2] * db[2];
  x = 2 / (length * x);
  for (int c = 0; c < 3; ++c) {
    const int i = permutation[c][0];
    const int j = permutation[c][1];
    const int k = permutation[c][2];
    // pa
    {
      Vec3 term1(real(0), real(0), real(0));
      Vec3 term2(real(0), real(0), real(0));
      Vec3 tmp(real(0), real(0), real(0));
      // first term
      dj::MulVecMatrix3x3<real>(db[k](), (real (*)[3]) &dajpi[j][0], term1());
      dj::MulVecMatrix3x3<real>(db[j](), (real (*)[3]) &dajpi[k][0], tmp());
      term1 = term1 - tmp;
      // second term
      for (int n = 0; n < 3; ++n) {
        dj::MulVecMatrix3x3<real>(db[n](), (real (*)[3]) &dajpi[n][0], tmp());
        term2 = term2 + tmp;
      }
      //      P(term1, term2, darboux_vector[i]);
      omega_pa[i] = (term1) - (real(0.5) * darboux_vector[i] * length) * (term2);
      omega_pa[i] *= (x * rod::bending_twising_stiffness[i]);
    }
    // pb
    {
      Vec3 term1(real(0), real(0), real(0));
      Vec3 term2(real(0), real(0), real(0));
      Vec3 tmp(real(0), real(0), real(0));
      // first term
      dj::MulVecMatrix3x3<real>(db[k](), (real (*)[3]) &dajpi[j][1], term1());
      dj::MulVecMatrix3x3<real>(db[j](), (real (*)[3]) &dajpi[k][1], tmp());
      term1 = term1 - tmp;
      // third term
      dj::MulVecMatrix3x3<real>(da[k](), (real (*)[3]) &dbjpi[j][0], tmp());
      term1 = term1 - tmp;
      dj::MulVecMatrix3x3<real>(da[j](), (real (*)[3]) &dbjpi[k][0], tmp());
      term1 = term1 + tmp;
      // second term
      for (int n = 0; n < 3; ++n) {
        dj::MulVecMatrix3x3<real>(db[n](), (real (*)[3]) &dajpi[n][1], tmp());
        term2 = term2 + tmp;
        dj::MulVecMatrix3x3<real>(da[n](), (real (*)[3]) &dbjpi[n][0], tmp());
        term2 = term2 + tmp;
      }
      omega_pb[i] = (term1) - (real(0.5) * darboux_vector[i] * length) * (term2);
      omega_pb[i] *= (x * rod::bending_twising_stiffness[i]);
    }
    // pc
    {
      Vec3 tmp(real(0), real(0), real(0));
      // first term
      Vec3 term1(real(0), real(0), real(0));
      dj::MulVecMatrix3x3<real>(da[k](), (real (*)[3]) &dbjpi[j][1], term1());
      dj::MulVecMatrix3x3<real>(da[j](), (real (*)[3]) &dbjpi[k][1], tmp());
      term1 = term1 - tmp;
      // second term
      Vec3 term2(real(0), real(0), real(0));
      for (int n = 0; n < 3; ++n) {
        dj::MulVecMatrix3x3<real>(da[n](), (real (*)[3]) &dbjpi[n][1], tmp());
        term2 = term2 + tmp;
      }
      omega_pc[i] = (term1) + (real(0.5) * darboux_vector[i] * length) * (term2);
      omega_pc[i] *= (-x * rod::bending_twising_stiffness[i]);
    }
    // pd
    {
      Vec3 tmp(real(0), real(0), real(0));
      // first term
      Vec3 term1(real(0), real(0), real(0));
      dj::MulVecMatrix3x3<real>(db[k](), (real (*)[3]) &dajpi[j][2], term1());
      dj::MulVecMatrix3x3<real>(db[j](), (real (*)[3]) &dajpi[k][2], tmp());
      term1 = term1 - tmp;
      // second term
      Vec3 term2(real(0), real(0), real(0));
      for (int n = 0; n < 3; ++n) {
        dj::MulVecMatrix3x3<real>(db[n](), (real (*)[3]) &dajpi[n][2], tmp());
        term2 = term2 + tmp;
      }
      omega_pd[i] = (term1) - (real(0.5) * darboux_vector[i] * length) * (term2);
      omega_pd[i] *= (x * rod::bending_twising_stiffness[i]);
    }
    // pe
    {
      Vec3 tmp(real(0), real(0), real(0));
      // first term
      Vec3 term1(real(0), real(0), real(0));
      dj::MulVecMatrix3x3<real>(da[k](), (real (*)[3]) &dbjpi[j][2], term1());
      //      if (i == 2) P(term1);
      dj::MulVecMatrix3x3<real>(da[j](), (real (*)[3]) &dbjpi[k][2], tmp());
      //      if (i == 2) P(tmp);
      term1 -= tmp;
      // second term
      Vec3 term2(real(0), real(0), real(0));
      for (int n = 0; n < 3; ++n) {
        //        if (i == 2) P(tmp);
        dj::MulVecMatrix3x3<real>(da[n](), (real (*)[3]) &dbjpi[n][2][0][0], tmp());
        //        if (i == 2) {
        //          P(da[n]);
        //          P(dbjpi[n][2]);
        //          P(tmp);
        //        }
        term2 += tmp;
      }
      //      if (i == 2) {
      //        P(i, j, k);
      //        P(term1, term2);
      //        P(darboux_vector[i], x);
      //      }
      omega_pe[i] = (term1) + (real(0.5) * darboux_vector[i] * length) * (term2);
      //      P((real(0.5) * darboux_vector[i] * avg_edge_length_) * (term2));
      //      P(omega_pe[i]);
      omega_pe[i] *= (-x * rod::bending_twising_stiffness[i]);
    }
  }

}

inline void PBDElasticRod::CorrectGravityForGhostPoints(real dt)
{
  OMP_FOR
  for (int e = 0; e < e_num_; ++e) {
    vm_[e] = real(0.5) * (vel_vert_[e] + vel_vert_[e + 1]);
  }
  real inv_dt = 1 / dt;
  Vec3 dt_g = dt * gravity_;
  Vec3 g_over_g_mag_2 = (1 / gravity_.MagnitudeSquare()) * gravity_;
  for (int e = 0; e < e_num_; ++e) {
    Vec3 accel = (vm_[e] - prev_vm_[e]) * inv_dt;
    real r = accel * g_over_g_mag_2;
    Vec3 correction = (1 - r) * dt_g;
    vel_edge_[e] -= correction;
    correction *= real(0.5);
    vel_vert_[e] += correction;
    vel_vert_[e + 1] += correction;
  }
  prev_vm_.swap(vm_);
}

void PBDElasticRod::FixEdge(int e)
{
  vert_[e] = rest_shape_.vert_[e];
  vert_[e + 1] = rest_shape_.vert_[e + 1];
  edge_vert_[e] = rest_shape_.edge_vert_[e];
}

void PBDElasticRod::ComputeMaterialFrame(ElasticRodShape * rod)
{
  for (unsigned e = 0; e < rod->edge_vert_.size(); ++e) {
    ComputeOneMaterialFrame(vert_[e], vert_[e + 1], edge_vert_[e], d_[e]);
  }
}

void PBDElasticRod::Render()
{
  super::Render();
}
